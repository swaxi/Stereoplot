# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Stereoplot
 Interactive stereographic projection plotting and analysis tool for QGIS
                             -------------------
        begin               : 2026-06-02
        copyright           : (C) 2026 by Julien Perret and Mark Jessell
        email               : julien.perret@uwa.edu.au; mark.jessell@uwa.edu.au
        git sha             : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                          *
 *   This program is free software; you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation; either version 2 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""

from qgis.PyQt.QtGui import *
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import to_hex, to_rgba
from matplotlib.path import Path
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.widgets import CheckButtons, Button
from collections import defaultdict
from .mplstereonet import *
from .mplstereonet import stereonet_math
from qgis.core import *
from qgis.gui import *
import os
from qgis.core import QgsProject
from math import asin,sin,degrees,radians,cos,tan,atan
import json
import re

# Use Arial for all Matplotlib-rendered plot text when available.
# Matplotlib will fall back to DejaVu Sans if Arial is not installed.
plt.rcParams.update({
    'font.family': 'Arial',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
})

class _BoundedLassoSelector:
    """Lasso selector that works entirely in axes coordinates (0-1 space).

    This avoids the numerically unstable inverse Lambert projection transform.
    The stereonet boundary is always the circle centred at (0.5, 0.5) with
    radius 0.5 in axes coordinates.  When the cursor leaves that circle the
    path is clipped to the boundary; arc vertices are interpolated so the line
    visually hugs the edge.  The onselect callback receives a list of
    (ax_x, ax_y) vertices in axes coordinates.
    """

    _R = 0.5  # stereonet radius in axes coords

    def __init__(self, ax, onselect):
        self.ax = ax
        self.onselect = onselect
        self._active = False
        self._verts = []          # list of (ax_x, ax_y)
        self._last_clipped = False
        self._last_angle = 0.0
        # Line drawn in axes-coordinate space so no data transform is needed
        self._line, = ax.plot([], [], color='black', linewidth=0.8,
                              transform=ax.transAxes)
        self._line.set_visible(False)
        canvas = ax.figure.canvas
        self._cids = [
            canvas.mpl_connect('button_press_event', self._on_press),
            canvas.mpl_connect('button_release_event', self._on_release),
            canvas.mpl_connect('motion_notify_event', self._on_move),
        ]

    def _disp_to_axes(self, x, y):
        """Convert a single display-coord point to axes coords."""
        return self.ax.transAxes.inverted().transform([[x, y]])[0]

    def _clip(self, x, y):
        """Convert display coords to axes coords, clipped to the stereonet circle.
        Returns (ax_x, ax_y, is_clipped, angle_from_centre)."""
        ax_x, ax_y = self._disp_to_axes(x, y)
        dx, dy = ax_x - 0.5, ax_y - 0.5
        angle = np.arctan2(dy, dx)
        if np.hypot(dx, dy) <= self._R:
            return ax_x, ax_y, False, angle
        return (0.5 + self._R * np.cos(angle),
                0.5 + self._R * np.sin(angle),
                True, angle)

    def _arc_verts(self, from_angle, to_angle):
        """Axes-coord vertices on the boundary arc, start-exclusive, end-inclusive."""
        diff = (to_angle - from_angle + np.pi) % (2 * np.pi) - np.pi
        n = max(2, int(abs(diff) * 20))
        angles = np.linspace(from_angle, from_angle + diff, n + 1)[1:]
        return [(0.5 + self._R * np.cos(a), 0.5 + self._R * np.sin(a))
                for a in angles]

    def _update_line(self):
        if self._verts:
            xs, ys = zip(*self._verts)
        else:
            xs, ys = [], []
        self._line.set_data(xs, ys)

    def _on_press(self, event):
        if event.button != 1:
            return
        bbox = self.ax.get_window_extent()
        if not (bbox.x0 <= event.x <= bbox.x1 and bbox.y0 <= event.y <= bbox.y1):
            return
        self._active = True
        self._verts = []
        ax_x, ax_y, clipped, angle = self._clip(event.x, event.y)
        self._verts.append((ax_x, ax_y))
        self._last_clipped = clipped
        self._last_angle = angle
        self._line.set_visible(True)
        self._update_line()
        self.ax.figure.canvas.draw_idle()

    def _on_move(self, event):
        if not self._active:
            return
        ax_x, ax_y, clipped, angle = self._clip(event.x, event.y)
        if self._last_clipped and clipped:
            self._verts.extend(self._arc_verts(self._last_angle, angle))
            self._last_angle = angle
        else:
            self._verts.append((ax_x, ax_y))
            self._last_clipped = clipped
            if clipped:
                self._last_angle = angle
        self._update_line()
        self.ax.figure.canvas.draw_idle()

    def _on_release(self, event):
        if not self._active or event.button != 1:
            return
        self._active = False
        self._line.set_visible(False)
        self.ax.figure.canvas.draw_idle()
        self.onselect(self._verts)

    def disconnect(self):
        for cid in self._cids:
            self.ax.figure.canvas.mpl_disconnect(cid)


class StereonetSettingsDialog(QDialog):
    """Settings dialog for controlling stereonet plot style.

    If config_path points to an existing stereonet.json the dialog reads from
    and writes back to that file.  Otherwise values are persisted in QSettings.
    """

    _QSETTINGS_ORG = 'Stereoplot'
    _QSETTINGS_APP = 'Stereoplot'
    _DEFAULTS = {
        'showGtCircles': False, 'showContours': True,
        'showKinematics': False, 'linPlanes': True, 'roseDiagram': False,
        'fitGirdle': False, 'dataType': 'Planes Only',
        'kinematicsField': None,
        'kinematicsAnchor': 'Plane pole',
        'classificationEnabled': False,
        'classificationField': None,
        'filterEnabled': False,
        'filterExpression': '',
        'classificationLayerSignature': '',
        'plotLayerSignature': '',
    }

    def __init__(self, parent=None, config_path=None, detected_data_type=None,
                 selected_layers=None, kinematics_candidate_fields=None):
        super().__init__(parent)
        self._config_path = config_path
        self._selected_layers = selected_layers or []
        self._kinematics_candidate_fields = kinematics_candidate_fields or []
        self._selected_kinematics_field = None
        self._filter_expression = ''
        self.setWindowTitle('Stereographic Projection Settings')
        self.setModal(True)

        cfg = self._load()

        current_layer_signature = self._selected_layer_signature()
        saved_plot_layer_signature = str(
            cfg.get('plotLayerSignature') or cfg.get('classificationLayerSignature', '') or '')
        # Rose Diagram is a plot-mode choice, not an intrinsic property of a
        # layer.  Do not carry it across to a newly selected layer, otherwise
        # all stereonet-specific tools remain disabled even when the layer has
        # been correctly detected as planar/linear/combined.
        if (cfg.get('roseDiagram', False) and detected_data_type and
                current_layer_signature and
                current_layer_signature != saved_plot_layer_signature):
            cfg['roseDiagram'] = False

        outer = QVBoxLayout()
        outer.addWidget(QLabel(
            'Select Features to Plot '
            '(Lower Hemisphere, Equal-Area Stereonet Projection):'))

        row = QHBoxLayout()
        self.gtCircles_cb  = QCheckBox('Great Circles')
        self.contours_cb   = QCheckBox('Contours')
        self.linPlanes_cb  = QCheckBox('Lineation-bearing Planes')
        self.rose_cb       = QCheckBox('Rose Diagram')
        self.kinematics_cb = QCheckBox('Kinematics')
        self.fitGirdle_cb  = QCheckBox('Best Fit Girdle')

        self.gtCircles_cb.setChecked( cfg['showGtCircles'])
        self.contours_cb.setChecked(  cfg['showContours'])
        self.linPlanes_cb.setChecked( cfg['linPlanes'])
        self.rose_cb.setChecked(      cfg['roseDiagram'])
        self.fitGirdle_cb.setChecked( cfg['fitGirdle'])

        self._selected_kinematics_field = cfg.get('kinematicsField')
        kin_available = self._kinematics_context_available()
        self.kinematics_cb.setEnabled(kin_available)
        self.kinematics_cb.setChecked(bool(cfg.get('showKinematics', False)) and kin_available)
        self.kinematics_cb.toggled.connect(self._on_kinematics_toggled)

        for cb in [self.gtCircles_cb, self.contours_cb, self.linPlanes_cb,
                   self.rose_cb, self.kinematics_cb, self.fitGirdle_cb]:
            row.addWidget(cb)

        update_btn = QPushButton('Update Settings')
        update_btn.clicked.connect(self._save_and_close)
        row.addStretch()
        row.addWidget(update_btn)

        outer.addLayout(row)

        dt_row = QHBoxLayout()
        dt_row.addWidget(QLabel('Data to plot:'))
        self.dataType_cb = QComboBox()
        self.dataType_cb.addItems([
            'Planes Only',
            'Lineations Only',
            'Lineations with Bearing Planes'
        ])

        initial_data_type = detected_data_type or cfg.get('dataType', 'Planes Only')
        if initial_data_type == 'Lineations with Planes':
            initial_data_type = 'Lineations with Bearing Planes'
        self.dataType_cb.setCurrentText(initial_data_type)

        self.dataType_cb.currentTextChanged.connect(self._on_data_type_changed)
        dt_row.addWidget(self.dataType_cb)
        dt_row.addStretch()
        outer.addLayout(dt_row)

        kin_anchor_row = QHBoxLayout()
        kin_anchor_row.addWidget(QLabel('Position of hangingwall displacement arrow:'))
        self.kinematics_anchor_cb = QComboBox()
        self.kinematics_anchor_cb.addItems(['Plane pole', 'Lineation'])
        saved_anchor = cfg.get('kinematicsAnchor', 'Plane pole')
        if saved_anchor not in ('Plane pole', 'Lineation'):
            saved_anchor = 'Plane pole'
        self.kinematics_anchor_cb.setCurrentText(saved_anchor)
        self.kinematics_anchor_cb.setEnabled(self.kinematics_cb.isChecked())
        self.kinematics_cb.toggled.connect(self.kinematics_anchor_cb.setEnabled)
        kin_anchor_row.addWidget(self.kinematics_anchor_cb)
        kin_anchor_row.addStretch()
        outer.addLayout(kin_anchor_row)

        # Apply dependent tool availability only after all widgets referenced
        # by _on_data_type_changed() have been created.
        self._on_data_type_changed(self.dataType_cb.currentText())


        class_group = QGroupBox('Classification')
        class_layout = QVBoxLayout()
        self.classification_cb = QCheckBox('Enable classification')
        current_layer_signature = self._selected_layer_signature()
        saved_layer_signature = str(cfg.get('classificationLayerSignature', '') or '')
        # Preserve classification while the user changes the selected subset
        # within the same layer.  Reset it only when the selected layer changes.
        self.classification_cb.setChecked(
            bool(cfg.get('classificationEnabled', False)) and
            bool(current_layer_signature) and
            current_layer_signature == saved_layer_signature
        )
        class_layout.addWidget(self.classification_cb)
        class_field_row = QHBoxLayout()
        class_field_row.addWidget(QLabel('Classify by:'))
        self.classification_field_cb = QComboBox()
        self.classification_field_cb.addItem('', '')
        for field_name, label in self._available_attribute_field_items():
            self.classification_field_cb.addItem(label, field_name)
        saved_class_field = cfg.get('classificationField') or ''
        if saved_class_field and self._combo_index_by_data(self.classification_field_cb, saved_class_field) == -1:
            self.classification_field_cb.addItem(saved_class_field, saved_class_field)
        saved_index = self._combo_index_by_data(self.classification_field_cb, saved_class_field)
        if saved_index != -1:
            self.classification_field_cb.setCurrentIndex(saved_index)
        self.classification_field_cb.setEnabled(self.classification_cb.isChecked())
        self.classification_cb.toggled.connect(self.classification_field_cb.setEnabled)
        class_field_row.addWidget(self.classification_field_cb)
        class_field_row.addStretch()
        class_layout.addLayout(class_field_row)
        class_group.setLayout(class_layout)
        outer.addWidget(class_group)

        filter_group = QGroupBox('Data Filter')
        filter_layout = QVBoxLayout()
        self.filter_cb = QCheckBox('Enable filtering')
        self.filter_cb.setChecked(bool(cfg.get('filterEnabled', False)))
        filter_layout.addWidget(self.filter_cb)
        filter_row = QHBoxLayout()
        self.filter_expr_le = QLineEdit(cfg.get('filterExpression', '') or '')
        self.filter_expr_le.setPlaceholderText("Example: \"Generation\" = '1' AND \"Kinematics\" IS NOT NULL")
        self.filter_expr_le.setEnabled(self.filter_cb.isChecked())
        filter_row.addWidget(self.filter_expr_le)
        self.filter_build_btn = QPushButton('Build Filter...')
        self.filter_build_btn.setEnabled(self.filter_cb.isChecked())
        self.filter_build_btn.clicked.connect(self._build_filter_expression)
        filter_row.addWidget(self.filter_build_btn)
        self.filter_cb.toggled.connect(self.filter_expr_le.setEnabled)
        self.filter_cb.toggled.connect(self.filter_build_btn.setEnabled)
        filter_layout.addLayout(filter_row)
        filter_group.setLayout(filter_layout)
        outer.addWidget(filter_group)

        self._rose_exclusive_widgets = [
            self.gtCircles_cb,
            self.contours_cb,
            self.linPlanes_cb,
            self.kinematics_cb,
            self.fitGirdle_cb,
        ]
        self.rose_cb.toggled.connect(self._on_rose_diagram_toggled)
        self._on_rose_diagram_toggled(self.rose_cb.isChecked())

        self.setLayout(outer)


    def _on_rose_diagram_toggled(self, checked):
        """Make Rose Diagram mutually exclusive with stereonet-specific options.

        Rose diagrams use azimuthal frequencies only. Great circles, contours,
        lineation-bearing planes, kinematic arrows and best-fit girdles are
        stereonet-specific overlays, so they are disabled while the rose
        diagram mode is active.
        """
        if checked:
            for widget in self._rose_exclusive_widgets:
                widget.blockSignals(True)
                widget.setChecked(False)
                widget.setEnabled(False)
                widget.blockSignals(False)
            self.kinematics_anchor_cb.setEnabled(False)
            return

        # Restore normal availability rules when returning to stereonet mode.
        self.gtCircles_cb.setEnabled(True)
        self.contours_cb.setEnabled(True)
        self.fitGirdle_cb.setEnabled(True)
        kin_available = self._kinematics_context_available()
        self.kinematics_cb.setEnabled(kin_available)
        self.kinematics_anchor_cb.setEnabled(self.kinematics_cb.isChecked() and kin_available)
        self._on_data_type_changed(self.dataType_cb.currentText())


    def _available_attribute_field_items(self):
        """Return unique attribute fields as (field_name, display_label).

        The display label includes the QGIS field alias when available, while
        the stored combo-box value remains the real provider field name.
        """
        items = []
        seen = set()
        for layer in self._selected_layers:
            if layer.type() != QgsMapLayer.VectorLayer:
                continue
            for index, field in enumerate(layer.fields()):
                name = field.name()
                if name in seen:
                    continue
                seen.add(name)
                alias = ''
                try:
                    alias = layer.attributeAlias(index) or ''
                except Exception:
                    alias = ''
                label = alias if alias else name
                items.append((name, label))
        return items

    def _available_attribute_fields(self):
        """Return all unique real attribute field names from selected vector layers."""
        return [field_name for field_name, _ in self._available_attribute_field_items()]

    def _selected_layer_signature(self):
        """Return a stable signature for the selected vector layer set.

        Classification is preserved while working on the same layer, even if
        the selected features change.  It is reset only when the user switches
        to a different selected vector layer or layer combination.
        """
        signatures = []
        for layer in self._selected_layers:
            if layer.type() != QgsMapLayer.VectorLayer:
                continue
            try:
                signatures.append(layer.id())
            except Exception:
                signatures.append(layer.name())
        return '|'.join(sorted(str(item) for item in signatures))

    @staticmethod
    def _combo_index_by_data(combo, value):
        for index in range(combo.count()):
            if combo.itemData(index) == value:
                return index
        return -1

    @staticmethod
    def _combo_current_data(combo):
        value = combo.currentData()
        return value if value is not None else combo.currentText().strip()

    def _expression_layer(self):
        """Return the first selected vector layer for the expression builder."""
        for layer in self._selected_layers:
            if layer.type() == QgsMapLayer.VectorLayer:
                return layer
        return None

    def _build_filter_expression(self):
        """Open the native QGIS layer-filter dialog where available.

        QgsQueryBuilder is the same family of dialog used by QGIS for layer
        filtering.  If it is unavailable in a given QGIS build, fall back to the
        expression builder so the plugin remains usable.
        """
        layer = self._expression_layer()
        if layer is None:
            QMessageBox.warning(self, 'Data Filter', 'No selected vector layer is available for filter construction.')
            return
        current_expression = self.filter_expr_le.text().strip()

        try:
            dlg = QgsQueryBuilder(layer, self)
            if current_expression and hasattr(dlg, 'setSql'):
                dlg.setSql(current_expression)
            if dlg.exec():
                if hasattr(dlg, 'sql'):
                    expression = dlg.sql()
                elif hasattr(dlg, 'sqlText'):
                    expression = dlg.sqlText()
                else:
                    expression = current_expression
                self.filter_expr_le.setText(expression)
            return
        except Exception:
            pass

        try:
            dlg = QgsExpressionBuilderDialog(layer, current_expression, self)
        except TypeError:
            dlg = QgsExpressionBuilderDialog(layer, current_expression, self, 'generic')
        if dlg.exec():
            expression = dlg.expressionText()
            self.filter_expr_le.setText(expression)

    def _validate_filter_expression(self):
        if not self.filter_cb.isChecked():
            return True
        expression_text = self.filter_expr_le.text().strip()
        if not expression_text:
            QMessageBox.critical(self, 'Data Filter Error', 'Filtering is enabled, but no filter expression has been provided.')
            return False
        expression = QgsExpression(expression_text)
        if expression.hasParserError():
            QMessageBox.critical(self, 'Data Filter Error', expression.parserErrorString())
            return False
        return True

    def _validate_classification_request(self):
        if not self.classification_cb.isChecked():
            return True
        field_name = self._combo_current_data(self.classification_field_cb).strip()
        if not field_name:
            QMessageBox.critical(self, 'Classification Error', 'Classification is enabled, but no classification field has been selected.')
            return False
        if field_name not in self._available_attribute_fields():
            QMessageBox.critical(self, 'Classification Error', 'The selected classification field does not exist in the selected layer(s).')
            return False
        return True

    def _on_data_type_changed(self, text):
        """Synchronise tool availability with the selected plotting mode.

        The user can deliberately override the automatically detected data
        type.  The dependent options must therefore be refreshed immediately
        from the dropdown value, not only when the settings dialog is reopened.
        """
        if getattr(self, 'rose_cb', None) is not None and self.rose_cb.isChecked():
            return

        is_planes = text == 'Planes Only'
        is_lines = text == 'Lineations Only'
        is_combined = text == 'Lineations with Bearing Planes'

        if is_combined:
            self.linPlanes_cb.setChecked(True)
            self.linPlanes_cb.setEnabled(True)
        else:
            self.linPlanes_cb.setChecked(False)
            self.linPlanes_cb.setEnabled(False)

        self.gtCircles_cb.setEnabled(is_planes or is_combined)
        self.contours_cb.setEnabled(True)
        self.fitGirdle_cb.setEnabled(True)

        kin_available = self._kinematics_context_available() and (is_lines or is_combined)
        self.kinematics_cb.setEnabled(kin_available)
        anchor_widget = getattr(self, 'kinematics_anchor_cb', None)
        if not kin_available:
            self.kinematics_cb.blockSignals(True)
            self.kinematics_cb.setChecked(False)
            self.kinematics_cb.blockSignals(False)
            if anchor_widget is not None:
                anchor_widget.setEnabled(False)
        elif anchor_widget is not None:
            anchor_widget.setEnabled(self.kinematics_cb.isChecked())

    @staticmethod
    def _normalise_token(value):
        return re.sub(r'[^a-z0-9]+', '', str(value or '').lower())

    @classmethod
    def _normalise_kinematic_value(cls, value):
        lookup = {
            'strikeslip': 'strike-slip',
            'strike': 'strike-slip',
            'sinslip': 'sinistral',
            'sinistral': 'sinistral',
            'sinistralslip': 'sinistral',
            'leftlateral': 'sinistral',
            'leftlateralslip': 'sinistral',
            'leftlatslip': 'sinistral',
            'leftlat': 'sinistral',
            'sin': 'sinistral',
            'dextral': 'dextral',
            'dextralslip': 'dextral',
            'rightlateral': 'dextral',
            'rightlateralslip': 'dextral',
            'rightlatslip': 'dextral',
            'rightlat': 'dextral',
            'dex': 'dextral',
            'dipslip': 'dip-slip',
            'normal': 'normal',
            'normalslip': 'normal',
            'normalfault': 'normal',
            'normalfaulting': 'normal',
            'extensional': 'normal',
            'extension': 'normal',
            'reverse': 'reverse',
            'reverseslip': 'reverse',
            'reversefault': 'reverse',
            'reversefaulting': 'reverse',
            'thrust': 'reverse',
            'thrustslip': 'reverse',
            'thrustfault': 'reverse',
            'compressional': 'reverse',
            'compression': 'reverse',
        }
        return lookup.get(cls._normalise_token(value))

    @staticmethod
    def _field_exists_on_layer(layer, field_name):
        return layer.fields().lookupField(field_name) != -1

    def _layer_has_lineation(self, layer):
        az_names = ['Azimuth', 'azimuth', 'Bearing', 'bearing', 'Trend', 'TREND']
        pl_names = ['Plunge', 'plunge']
        has_az = any(self._field_exists_on_layer(layer, name) for name in az_names)
        has_pl = any(self._field_exists_on_layer(layer, name) for name in pl_names)
        return has_az and has_pl

    def _layer_has_bearing_plane(self, layer):
        strike_names = ['Strike_RHR', 'Strike', 'strike', 'Strike_ref', 'Strike_Ref', 'strike_ref']
        dipdir_names = ['Dip_Direction', 'Dip_Dir', 'DipDirection', 'dip_direction',
                        'DipDir', 'DIPDIR', 'DipDir_ref', 'Dip_Dir_ref',
                        'DipDirection_ref', 'Dip_Direction_ref', 'DipDirection_Ref',
                        'dipdir_ref', 'dip_dir_ref', 'dip_direction_ref']
        dip_names = ['Dip', 'dip', 'Dip_ref', 'Dip_Ref', 'dip_ref']
        has_orientation = any(self._field_exists_on_layer(layer, name) for name in strike_names + dipdir_names)
        has_dip = any(self._field_exists_on_layer(layer, name) for name in dip_names)
        return has_orientation and has_dip

    def _kinematics_context_available(self):
        has_lineation = any(self._layer_has_lineation(layer) for layer in self._selected_layers
                            if layer.type() == QgsMapLayer.VectorLayer)
        return has_lineation and bool(self._kinematics_candidate_fields)

    def _on_kinematics_toggled(self, checked):
        if not checked:
            return
        if self._selected_kinematics_field in self._kinematics_candidate_fields:
            return
        if not self._kinematics_candidate_fields:
            QMessageBox.warning(self, 'Kinematics',
                                'No valid kinematics attribute field could be identified.')
            self.kinematics_cb.blockSignals(True)
            self.kinematics_cb.setChecked(False)
            self.kinematics_cb.blockSignals(False)
            return
        field, ok = QInputDialog.getItem(
            self, 'Kinematics Field',
            'Select the attribute field containing kinematic sense values:',
            self._kinematics_candidate_fields, 0, False)
        if ok and field:
            self._selected_kinematics_field = field
        else:
            self.kinematics_cb.blockSignals(True)
            self.kinematics_cb.setChecked(False)
            self.kinematics_cb.blockSignals(False)

    def _field_has_recognised_kinematic_values(self, field_name):
        """Return True if the field contains recognised kinematic values.

        Selected features are checked first, but the whole layer is inspected as
        a fallback.  This avoids disabling kinematic plotting when the current
        selection only contains non-plottable values such as ``Unknown`` or when
        QGIS has not yet propagated the selection state to the dialog.
        """
        for layer in self._selected_layers:
            if layer.type() != QgsMapLayer.VectorLayer:
                continue
            if not self._field_exists_on_layer(layer, field_name):
                continue

            features = list(layer.selectedFeatures())
            if not features:
                try:
                    features = list(layer.getFeatures())
                except Exception:
                    features = []

            for feature in features:
                value = _attr(feature[field_name])
                if value is not None and self._normalise_kinematic_value(value):
                    return True
        return False

    def _validate_kinematics_request(self):
        if self._selected_kinematics_field not in self._kinematics_candidate_fields:
            QMessageBox.critical(self, 'Kinematics Error',
                                 'No valid kinematics attribute field could be identified.')
            return False
        if not any(self._field_exists_on_layer(layer, self._selected_kinematics_field)
                   for layer in self._selected_layers
                   if layer.type() == QgsMapLayer.VectorLayer):
            QMessageBox.critical(self, 'Kinematics Error',
                                 'No valid kinematics attribute field could be identified.')
            return False
        if not self._field_has_recognised_kinematic_values(self._selected_kinematics_field):
            QMessageBox.critical(
                self, 'Kinematics Error',
                'The selected kinematics field does not contain recognised kinematic values '
                '(Sinistral, Dextral, Normal, Reverse, etc.).')
            return False
        # Feature-level bearing-plane values are checked during plotting:
        # kinematic arrows are generated only where a valid associated plane
        # orientation is available.  Features with NULL strike/dip-direction/dip
        # values are treated as non-plottable for kinematic arrows.
        return True

    def _load(self):
        """Load config from JSON file if present, else QSettings, else defaults."""
        if self._config_path and os.path.exists(self._config_path):
            with open(self._config_path, 'r') as f:
                cfg = json.load(f)
            return {k: cfg.get(k, v) for k, v in self._DEFAULTS.items()}
        s = QSettings(self._QSETTINGS_ORG, self._QSETTINGS_APP)
        if s.contains('showContours'):
            result = {}
            for k, v in self._DEFAULTS.items():
                result[k] = s.value(k, v, type=bool) if isinstance(v, bool) else s.value(k, v)
            return result
        return dict(self._DEFAULTS)

    def _save_and_close(self):
        rose_mode = self.rose_cb.isChecked()
        if rose_mode:
            self.gtCircles_cb.setChecked(False)
            self.contours_cb.setChecked(False)
            self.linPlanes_cb.setChecked(False)
            self.kinematics_cb.setChecked(False)
            self.fitGirdle_cb.setChecked(False)
        elif self.dataType_cb.currentText() == 'Lineations with Bearing Planes':
            self.linPlanes_cb.setChecked(True)
        else:
            self.linPlanes_cb.setChecked(False)

        if (not rose_mode and self.kinematics_cb.isChecked() and
                not self._validate_kinematics_request()):
            return
        if not self._validate_classification_request():
            return
        if not self._validate_filter_expression():
            return

        cfg = {
            'showGtCircles':  False if rose_mode else self.gtCircles_cb.isChecked(),
            'showContours':   False if rose_mode else self.contours_cb.isChecked(),
            'showKinematics': False if rose_mode else self.kinematics_cb.isChecked(),
            'linPlanes':      False if rose_mode else self.linPlanes_cb.isChecked(),
            'roseDiagram':    rose_mode,
            'fitGirdle':      False if rose_mode else self.fitGirdle_cb.isChecked(),
            'dataType':       self.dataType_cb.currentText(),
            'kinematicsField': self._selected_kinematics_field,
            'kinematicsAnchor': self.kinematics_anchor_cb.currentText(),
            'classificationEnabled': self.classification_cb.isChecked(),
            'classificationField': self._combo_current_data(self.classification_field_cb).strip() or None,
            'filterEnabled': self.filter_cb.isChecked(),
            'filterExpression': self.filter_expr_le.text().strip(),
            'classificationLayerSignature': self._selected_layer_signature(),
            'plotLayerSignature': self._selected_layer_signature(),
        }
        if self._config_path:
            os.makedirs(os.path.dirname(self._config_path), exist_ok=True)
            with open(self._config_path, 'w') as f:
                json.dump(cfg, f, indent=4)
        else:
            s = QSettings(self._QSETTINGS_ORG, self._QSETTINGS_APP)
            for k, v in cfg.items():
                s.setValue(k, v)
        self.accept()

    @staticmethod
    def load_qsettings():
        """Return a stereoConfig dict from QSettings, or None if not yet saved."""
        s = QSettings(StereonetSettingsDialog._QSETTINGS_ORG,
                      StereonetSettingsDialog._QSETTINGS_APP)
        if not s.contains('showContours'):
            return None
        result = {}
        for k, v in StereonetSettingsDialog._DEFAULTS.items():
            result[k] = s.value(k, v, type=bool) if isinstance(v, bool) else s.value(k, v)
        return result


def _attr(v):
    """Return a QGIS feature attribute as a Python value, or None if null."""
    if v is None:
        return None
    if hasattr(v, 'isNull') and v.isNull():
        return None
    return v


# =============================================================================
# CATEGORY ORIENTATION STATISTICS OVERLAYS
# =============================================================================
# The interactive category panel (_open_category_panel) can reveal a hidden
# representative-orientation marker (and, optionally, a text label) per
# category, built once alongside the normal per-sample poles/lines, so
# toggling "Show mean or Kamb maximum orientation" is instant (no
# re-plotting) while the underlying points stay clickable. Two alternative
# statistics are built for each category, and the panel lets the user pick
# which one is shown:
#
#   'mean'        - the true spherical (Fisher) mean vector, not a naive
#                    arithmetic average of angles, via find_mean_vector().
#   'contour_max' - the peak of a per-category Kamb density grid (the same
#                    exponential-smoothing method/sigma as the plot's own
#                    density-contour overlay), via density_grid(). This can
#                    differ meaningfully from the mean for multi-modal or
#                    non-Fisherian scatter.

def _orientation_position(kind, mode, values1, values2):
    """Return (v1, v2): strike/dip (kind='planes') or plunge/bearing
    (kind='lines') of the requested representative orientation ('mean' or
    'contour_max') for one category's measurements."""
    measurement = 'poles' if kind == 'planes' else 'lines'
    if mode == 'contour_max':
        lon_grid, lat_grid, z = density_grid(
            values1, values2, measurement=measurement,
            method='exponential_kamb', sigma=1.5, gridsize=50)
        peak = np.unravel_index(np.argmax(z), z.shape)
        lon, lat = float(lon_grid[peak]), float(lat_grid[peak])
        if kind == 'planes':
            strike_arr, dip_arr = stereonet_math.geographic2pole(lon, lat)
            return float(strike_arr[0]), float(dip_arr[0])
        plunge_arr, bearing_arr = stereonet_math.geographic2plunge_bearing(lon, lat)
        return float(plunge_arr[0]), float(bearing_arr[0])

    mean_vector, _r_value = find_mean_vector(values1, values2, measurement=measurement)
    if kind == 'planes':
        strike_arr, dip_arr = plunge_bearing2pole(*mean_vector)
        return float(strike_arr[0]), float(dip_arr[0])
    return float(mean_vector[0]), float(mean_vector[1])


def _build_orientation_stat_artists(ax, kind, values1, values2, style):
    """Build hidden representative-orientation markers, text labels and (for
    'planes' only) great circles for one category's poles ('planes') or
    lines ('lines') - one set for each of the 'mean' and 'contour_max'
    statistics (see module docstring above).

    All artists are created hidden (visible=False); the interactive category
    panel reveals whichever statistic is selected and fades the raw points
    on a "Show mean or Kamb maximum orientation" toggle. Each marker's size
    follows the category's own markersize (scaled up slightly) and its
    colour/marker/fill are kept in sync with the category style whenever it
    changes; each label defaults to the same colour, with a white "mask"
    behind the text so it stays legible over other plotted lines/points. For
    planes, the pole marker is a shorthand for a plane orientation, so its
    label reads "Pole to <strike>/<dip>" and the plane itself is also drawn,
    as a dashed great circle through that pole.

    Returns (stats_by_mode, label_by_mode, circle_by_mode): all {'mean': ...,
    'contour_max': ...} dicts, or (None, None, None) if there is no data.
    stats_by_mode[mode] is the list of Line2D marker artists from
    ax.pole()/ax.line(); label_by_mode[mode] is the matching Annotation;
    circle_by_mode[mode] is the list of Line2D great-circle artists from
    ax.plane() for kind='planes', or None for kind='lines'.
    """
    if not values1:
        return None, None, None

    colour = style.get('color', '#000000')
    hollow = style.get('fill') == 'hollow'
    face = 'white' if hollow else colour
    markersize = max(6.0, float(style.get('markersize', 5)) * 1.6)
    linewidth = max(1.2, float(style.get('linewidth', 1.0)) * 1.3)

    def _make_artist(v1, v2):
        if kind == 'planes':
            marker_artist = ax.pole(
                [v1], [v2], linestyle='none', marker=style.get('marker', 'o'),
                markerfacecolor=face, markeredgecolor=colour, markeredgewidth=1.4,
                markersize=markersize, zorder=16, visible=False)
            circle_artist = ax.plane(
                v1, v2, color=colour, linewidth=linewidth, linestyle='--',
                zorder=15, visible=False)
            text = f"Pole to {v1:03.0f}/{v2:02.0f}"
            lon, lat = stereonet_math.pole(v1, v2)
        else:
            marker_artist = ax.line(
                [v1], [v2], linestyle='none', marker=style.get('marker', 'o'),
                markerfacecolor=face, markeredgecolor=colour, markeredgewidth=1.4,
                markersize=markersize, zorder=16, visible=False)
            circle_artist = None
            text = f"{v1:02.0f}→{v2:03.0f}"
            lon, lat = stereonet_math.line(v1, v2)
        label_artist = ax.annotate(
            text, xy=(float(np.atleast_1d(lon)[0]), float(np.atleast_1d(lat)[0])),
            xytext=(8, 8), textcoords='offset points', color=colour,
            fontsize=8, fontfamily='Arial', ha='left', va='bottom',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.85),
            zorder=17, visible=False)
        return marker_artist, label_artist, circle_artist

    stats_by_mode = {}
    label_by_mode = {}
    circle_by_mode = {}
    for mode in ('mean', 'contour_max'):
        v1, v2 = _orientation_position(kind, mode, values1, values2)
        marker_artist, label_artist, circle_artist = _make_artist(v1, v2)
        stats_by_mode[mode] = marker_artist
        label_by_mode[mode] = label_artist
        circle_by_mode[mode] = circle_artist

    return stats_by_mode, label_by_mode, circle_by_mode


def _restyle_mean_artist(stats_by_mode, style):
    """Recolour/re-marker every hidden representative-orientation Line2D
    list (for both the 'mean' and 'contour_max' modes) to match a category's
    current style, leaving each artist's computed position untouched."""
    if not stats_by_mode:
        return
    colour = style.get('color', '#000000')
    face = 'white' if style.get('fill') == 'hollow' else colour
    markersize = max(6.0, float(style.get('markersize', 5)) * 1.6)
    for artist in stats_by_mode.values():
        if not artist:
            continue
        for item in artist:
            if hasattr(item, 'set_marker'):
                item.set_marker(style.get('marker', 'o'))
            if hasattr(item, 'set_markerfacecolor'):
                item.set_markerfacecolor(face)
            if hasattr(item, 'set_markeredgecolor'):
                item.set_markeredgecolor(colour)
            if hasattr(item, 'set_markersize'):
                item.set_markersize(markersize)


def _restyle_mean_circle(circle_by_mode, style):
    """Recolour a hidden mean/Kamb-maximum great circle (for both modes) to
    match a category's current style, leaving its computed orientation
    untouched."""
    if not circle_by_mode:
        return
    colour = style.get('color', '#000000')
    linewidth = max(1.2, float(style.get('linewidth', 1.0)) * 1.3)
    for artist in circle_by_mode.values():
        if not artist:
            continue
        for item in artist:
            if hasattr(item, 'set_color'):
                item.set_color(colour)
            if hasattr(item, 'set_linewidth'):
                item.set_linewidth(linewidth)


def classFactory(iface):
    return Stereonet(iface)

class Stereonet:
    def __init__(self, iface):
        self.iface = iface

    def initGui(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.contourAction = QAction(QIcon(str(dir_path)+"/icon.png"), u'Stereonet', self.iface.mainWindow())
        self.contourAction.triggered.connect(self.contourPlot)
        self.iface.addToolBarIcon(self.contourAction)

        self.settingsAction = QAction(
            QgsApplication.getThemeIcon('mActionOptions.svg'),
            u'Stereonet Settings', self.iface.mainWindow())
        self.settingsAction.triggered.connect(self.showSettings)
        self.iface.addToolBarIcon(self.settingsAction)

    def unload(self):
        self.iface.removeToolBarIcon(self.contourAction)
        del self.contourAction
        self.iface.removeToolBarIcon(self.settingsAction)
        del self.settingsAction

    def _structural_field_names(self):
        """Return accepted field-name variants for structural data detection."""
        return {
            'strike': ['Strike_RHR', 'Strike', 'strike'],
            'dipdir': ['Dip_Direction', 'Dip_Dir', 'DipDirection',
                       'dip_direction', 'DipDir', 'DIPDIR'],
            'dip': ['Dip', 'dip'],
            'azimuth': ['Azimuth', 'azimuth', 'Bearing', 'bearing',
                        'Trend', 'TREND'],
            'plunge': ['Plunge', 'plunge'],
            'strike_ref': ['Strike_ref', 'Strike_Ref', 'strike_ref'],
            'dip_ref': ['Dip_ref', 'Dip_Ref', 'dip_ref'],
            'dipdir_ref': ['DipDir_ref', 'Dip_Dir_ref', 'DipDirection_ref',
                           'Dip_Direction_ref', 'DipDirection_Ref',
                           'dipdir_ref', 'dip_dir_ref', 'dip_direction_ref'],
            'kinematics': ['Kinematics', 'kinematics', 'Kinematic', 'kinematic', 'Kin', 'kin', 'Movement', 'movement', 'SlipSense', 'Slip_Sense', 'slip_sense', 'ShearSense', 'Shear_Sense', 'shear_sense', 'SenseOfMovement', 'Sense_of_Movement', 'sense_of_movement'],
            'pitch_rhr': ['Pitch_RHR', 'Pitch_rhr', 'Pitch_Rhr', 'Pitch',
                          'pitch_rhr', 'RHR_pitch', 'rhr_pitch', 'pitch'],
        }

    def _field_exists(self, layer, fieldnames):
        for fieldname in fieldnames:
            field_index = layer.fields().lookupField(fieldname)
            if field_index != -1:
                return True, fieldname
        return False, False

    @staticmethod
    def _normalise_token(value):
        return re.sub(r'[^a-z0-9]+', '', str(value or '').lower())

    def _candidate_kinematics_fields(self, layers):
        """Return likely kinematics fields in selected vector layers.

        Detection is based on both field names and field contents.  The latter
        catches legacy or abbreviated datasets where the field alias/name is not
        one of the preferred variants but the values clearly contain kinematic
        classes such as Sinistral-slip, Reverse-slip or Unknown.
        """
        names = self._structural_field_names()['kinematics']
        accepted = {self._normalise_token(name) for name in names}
        candidates = []
        for layer in layers:
            if layer.type() != QgsMapLayer.VectorLayer:
                continue

            for field in layer.fields():
                field_name = field.name()
                norm = self._normalise_token(field_name)
                name_match = (
                    norm in accepted or
                    'kinematic' in norm or
                    'movement' in norm or
                    'slipsense' in norm or
                    'senseofmovement' in norm or
                    'shearsense' in norm or
                    (('sense' in norm or 'slip' in norm or 'shear' in norm) and
                     ('kin' in norm or 'move' in norm or 'slip' in norm or 'shear' in norm))
                )

                value_match = False
                if not name_match:
                    # Inspect a small sample of selected features first, then
                    # fall back to layer features.  Stop as soon as a recognised
                    # kinematic value is encountered.
                    try:
                        features = list(layer.selectedFeatures())
                    except Exception:
                        features = []
                    if not features:
                        try:
                            features = list(layer.getFeatures())
                        except Exception:
                            features = []
                    for i, feature in enumerate(features):
                        if i >= 200:
                            break
                        value = _attr(feature[field_name])
                        if value is not None and self._normalise_kinematic_value(value):
                            value_match = True
                            break

                if (name_match or value_match) and field_name not in candidates:
                    candidates.append(field_name)
        return candidates

    @classmethod
    def _normalise_kinematic_value(cls, value):
        lookup = {
            'strikeslip': 'strike-slip',
            'sinistral': 'sinistral',
            'sinistralslip': 'sinistral',
            'leftlateral': 'sinistral',
            'sin': 'sinistral',
            'dextral': 'dextral',
            'dextralslip': 'dextral',
            'rightlateral': 'dextral',
            'dex': 'dextral',
            'dipslip': 'dip-slip',
            'normal': 'normal',
            'normalslip': 'normal',
            'extensional': 'normal',
            'reverse': 'reverse',
            'reverseslip': 'reverse',
            'thrust': 'reverse',
            'compressional': 'reverse',
        }
        return lookup.get(cls._normalise_token(value))

    @staticmethod
    def _as_scalar(value):
        value = np.asarray(value).ravel()
        if len(value) == 0:
            return np.nan
        return float(value[0])

    def _append_kinematic_arrow_record(self, records, strike, dip, plunge,
                                       bearing, sense_value, category='All'):
        # Unknown/undetermined kinematic values are intentionally treated as
        # NULL: they are allowed in the dataset but do not generate arrows.
        sense = self._normalise_kinematic_value(sense_value)
        if sense not in ('sinistral', 'dextral', 'normal', 'reverse'):
            return
        if None in (strike, dip, plunge, bearing):
            return
        records.append({
            'strike': float(strike),
            'dip': float(dip),
            'plunge': float(plunge),
            'bearing': float(bearing),
            'sense': sense,
            'category': category,
        })

    def _plot_kinematic_arrows(self, ax, records, anchor='Plane pole', color='black'):
        """Draw kinematic arrows on either the plane pole or the lineation.

        The arrow anchor is controlled by the settings dialog:
        - ``Plane pole`` keeps the arrow tail at the pole to the bearing plane.
        - ``Lineation`` keeps the arrow tail at the plotted lineation point.

        Arrow direction follows the local tangent to the great circle passing
        through the plane pole and associated lineation. Arrow length is set
        once in stereonet/data coordinates. At the initial figure size it is
        calibrated to appear approximately 0.75 cm long; when the figure window
        is resized, the arrow scales naturally with the stereonet.
        """
        initial_arrow_length_cm = 0.75
        anchor = anchor if anchor in ('Plane pole', 'Lineation') else 'Plane pole'

        def _great_circle_point(start, end, fraction):
            xyz0 = np.asarray(stereonet_math.sph2cart(start[0], start[1]), dtype=float)
            xyz1 = np.asarray(stereonet_math.sph2cart(end[0], end[1]), dtype=float)
            dot = float(np.clip(np.dot(xyz0, xyz1), -1.0, 1.0))
            omega = np.arccos(dot)
            if np.isclose(omega, 0.0):
                return None
            sin_omega = np.sin(omega)
            xyz = ((np.sin((1.0 - fraction) * omega) / sin_omega) * xyz0 +
                   (np.sin(fraction * omega) / sin_omega) * xyz1)
            lon, lat = stereonet_math.cart2sph(*xyz)
            return np.array([float(lon), float(lat)])

        def _unit_direction_at_anchor(pole, line, anchor_point, sense):
            """Return display-space tangent direction for the selected anchor."""
            if anchor == 'Lineation':
                # Direction from pole towards lineation, evaluated at the
                # lineation end of the same pole-lineation great circle.
                neighbour = _great_circle_point(pole, line, 0.97)
                if neighbour is None or np.any(np.isnan(neighbour)):
                    return None
                anchor_disp = ax.transData.transform(anchor_point)
                neighbour_disp = ax.transData.transform(neighbour)
                direction = anchor_disp - neighbour_disp
            else:
                # Direction from pole towards lineation, evaluated at the pole.
                neighbour = _great_circle_point(pole, line, 0.03)
                if neighbour is None or np.any(np.isnan(neighbour)):
                    return None
                anchor_disp = ax.transData.transform(anchor_point)
                neighbour_disp = ax.transData.transform(neighbour)
                direction = neighbour_disp - anchor_disp

            norm = np.hypot(direction[0], direction[1])
            if np.isclose(norm, 0.0):
                return None
            direction = direction / norm

            if sense == 'normal':
                # Movement is from pole towards lineation. If anchored on the
                # lineation, the arrow continues along the same local tangent.
                pass
            elif sense == 'reverse':
                # Movement is from lineation towards pole.
                direction = -direction
            elif sense == 'sinistral':
                # Right-directed strike-slip arrow in display space.
                if direction[0] < 0:
                    direction = -direction
            elif sense == 'dextral':
                # Left-directed strike-slip arrow in display space.
                if direction[0] > 0:
                    direction = -direction
            else:
                return None
            return direction

        # Ensure transforms are initialised before converting the requested
        # initial display length to a stereonet/data-coordinate offset.
        ax.figure.canvas.draw_idle()

        arrow_artists = []

        for rec in records:
            pole_lon, pole_lat = mplstereonet.pole(rec['strike'], rec['dip'])
            line_lon, line_lat = mplstereonet.line(rec['plunge'], rec['bearing'])

            pole = np.array([self._as_scalar(pole_lon), self._as_scalar(pole_lat)])
            line = np.array([self._as_scalar(line_lon), self._as_scalar(line_lat)])
            if np.any(np.isnan(pole)) or np.any(np.isnan(line)):
                continue
            if np.allclose(pole, line):
                continue

            anchor_point = line if anchor == 'Lineation' else pole
            direction = _unit_direction_at_anchor(pole, line, anchor_point, rec['sense'])
            if direction is None:
                continue

            anchor_disp = ax.transData.transform(anchor_point)
            initial_length_px = (initial_arrow_length_cm / 2.54) * ax.figure.dpi
            end_disp = anchor_disp + direction * initial_length_px
            end_data = ax.transData.inverted().transform(end_disp)

            arrow = FancyArrowPatch(
                posA=tuple(anchor_point), posB=tuple(end_data),
                arrowstyle='-|>', mutation_scale=10,
                linewidth=1.0, color=color,
                shrinkA=0, shrinkB=0,
                transform=ax.transData, zorder=6)
            ax.add_patch(arrow)
            arrow_artists.append(arrow)

        return arrow_artists

    def _detect_data_type_from_layers(self, layers):
        """Infer the plotting mode from fields available in selected layers."""
        names = self._structural_field_names()
        has_plane = False
        has_line = False
        has_ref_plane = False

        for layer in layers:
            if layer.type() != QgsMapLayer.VectorLayer:
                continue

            strike_ok, _ = self._field_exists(layer, names['strike'])
            dipdir_ok, _ = self._field_exists(layer, names['dipdir'])
            dip_ok, _ = self._field_exists(layer, names['dip'])
            az_ok, _ = self._field_exists(layer, names['azimuth'])
            plunge_ok, _ = self._field_exists(layer, names['plunge'])
            sref_ok, _ = self._field_exists(layer, names['strike_ref'])
            ddref_ok, _ = self._field_exists(layer, names['dipdir_ref'])
            dref_ok, _ = self._field_exists(layer, names['dip_ref'])

            has_plane = has_plane or ((strike_ok or dipdir_ok) and dip_ok)
            has_line = has_line or (az_ok and plunge_ok)
            has_ref_plane = has_ref_plane or ((sref_ok or ddref_ok) and dref_ok)

        # Combined mode applies to linear data carrying either explicit
        # reference-plane fields (e.g. Strike_ref/Dip_ref or DipDir_ref/Dip_ref)
        # or regular planar fields on the same layer (e.g. Folds_PT with
        # Azimuth/Plunge plus Strike_RHR or Dip_Dir/Dip).
        if has_line and (has_ref_plane or has_plane):
            return 'Lineations with Bearing Planes'
        if has_line:
            return 'Lineations Only'
        if has_plane:
            return 'Planes Only'
        return None

    def showSettings(self):
        project_file = QgsProject.instance().fileName()
        config_path = None
        if project_file:
            config_path = os.path.join(
                os.path.dirname(os.path.abspath(project_file)),
                "99_COMMAND_FILES_PLUGIN",
                "stereonet.json")

        layers = self.iface.layerTreeView().selectedLayers()
        detected_data_type = self._detect_data_type_from_layers(layers)
        kinematics_candidate_fields = self._candidate_kinematics_fields(layers)

        dlg = StereonetSettingsDialog(
            self.iface.mainWindow(),
            config_path=config_path,
            detected_data_type=detected_data_type,
            selected_layers=layers,
            kinematics_candidate_fields=kinematics_candidate_fields)
        dlg.exec()
    
    def waxi_tangent_lineation_plot(self,ax,strikes, dips,kinematics,rhr,azs):
        """Makes a tangent lineation plot for normal faults with the given strikes,
        dips, and rakes."""
        sos=['Sinistral-slip','Dextral-slip','Normal-slip','Reverse-slip']

        for j,sos_it in enumerate(sos):

            dp=list()
            sp=list()
            rh=list()
            az=list()
            for i in range(len(strikes)):
                if(kinematics[i]==sos_it):

                    dp.append(dips[i])
                    sp.append(strikes[i])
                    rh.append(rhr[i])
                    if(sos_it in ['Normal-slip','Reverse-slip']):
                        if(azs[i]>90):
                            azimuth=180-azs[i]
                        else:
                            azimuth=azs[i]
                        if(rhr[i]=='YES'):
                            az.append(360-azimuth)
                        else:
                            az.append(180-azimuth)
                    else:
                        if(azs[i]>180 ):
                            azimuth=90-azs[i]
                        else:
                            azimuth=azs[i]
                        if(rhr[i]=='YES'):
                            azimuth=180-azimuth
                        else:
                            azimuth=360-azimuth
                        if(rhr[i]=='YES' and azs[i]>90):
                            azimuth=azs[i]
                        az.append(azimuth)                        


            # Calculate the position of the rake of the lineations, but don't plot yet
            rake_x, rake_y = mplstereonet.rake(sp, dp, az)
            
            # Calculate the direction the arrows should point
            # These are all normal faults, so the arrows point away from the center
            # Because we're plotting at the pole location, however, we need to flip this
            # from what we plotted with the "ball of string" plot.
            mag = np.hypot(rake_x, rake_y)
            u, v = -rake_x / mag, -rake_y / mag

            # Calculate the position of the poles
            pole_x, pole_y = mplstereonet.pole(sp, dp)
           
            # Plot the arrows centered on the pole locations...
            if(sos_it=='Sinistral-slip'):
                arrows = ax.quiver(pole_x, pole_y, u, v,  width=1, headwidth=4, units='dots', color='r',
                                pivot='tail')
            elif(sos_it=='Dextral-slip'):
                arrows = ax.quiver(pole_x, pole_y, -u, -v,  width=1, headwidth=4, units='dots', color='g',
                                pivot='tail')
            elif(sos_it=='Normal-slip'):
                arrows = ax.quiver(pole_x, pole_y, u,-v,  width=1, headwidth=4, units='dots', color='b',
                                pivot='tail')
            elif(sos_it=='Reverse-slip'):
                arrows = ax.quiver(pole_x, pole_y, -u,v,  width=1, headwidth=4, units='dots', color='m',
                                pivot='tail')


            #return arrows
    



    def _build_filter_expression(self, expression_text):
        """Return a compiled QgsExpression or None when filtering is disabled."""
        expression_text = (expression_text or '').strip()
        if not expression_text:
            return None
        expression = QgsExpression(expression_text)
        if expression.hasParserError():
            self.iface.messageBar().pushMessage(
                'Invalid stereonet filter expression',
                expression.parserErrorString(), level=Qgis.Warning, duration=8)
            return None
        return expression

    def _feature_passes_filter(self, layer, feature, expression):
        """Evaluate a compiled QGIS expression against one feature."""
        if expression is None:
            return True
        context = QgsExpressionContext()
        try:
            context.appendScopes(QgsExpressionContextUtils.globalProjectLayerScopes(layer))
        except Exception:
            pass
        context.setFeature(feature)
        value = expression.evaluate(context)
        if expression.hasEvalError():
            self.iface.messageBar().pushMessage(
                'Stereonet filter evaluation error',
                expression.evalErrorString() + ' — filter ignored for this feature.',
                level=Qgis.Warning, duration=8)
            return True
        return bool(value)

    def _category_value(self, feature, field_name):
        """Return a displayable category value, ignoring NULL/empty values."""
        if not field_name:
            return 'All'
        value = _attr(feature[field_name])
        if value is None or str(value).strip() == '':
            return None
        return str(value)

    @staticmethod
    def _default_category_style(index):
        """Return a simple default style for category index."""
        palette = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                   'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
        markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', '<', '>']
        base_colour = palette[index % len(palette)]
        return {
            'color': base_colour,
            'linecolor': base_colour,
            'arrowcolor': base_colour,
            'marker': markers[index % len(markers)],
            'markersize': 5,
            'linewidth': 1.0,
            'alpha': 1.0,
            'fill': 'full',
        }

    def _style_templates_path(self):
        """Return the project-level JSON file used for classification style templates."""
        project_file = QgsProject.instance().fileName()
        if project_file:
            config_dir = os.path.join(
                os.path.dirname(os.path.abspath(project_file)),
                '99_COMMAND_FILES_PLUGIN')
            return os.path.join(config_dir, 'stereonet_styles.json')
        return os.path.join(os.path.expanduser('~'), 'stereonet_styles.json')

    def _load_style_templates(self):
        """Load saved classification style templates from JSON."""
        path = self._style_templates_path()
        if not path or not os.path.exists(path):
            return {'version': 1, 'templates': {}}
        try:
            with open(path, 'r') as f:
                data = json.load(f)
        except Exception:
            return {'version': 1, 'templates': {}}
        if not isinstance(data, dict):
            return {'version': 1, 'templates': {}}
        data.setdefault('version', 1)
        data.setdefault('templates', {})
        if not isinstance(data['templates'], dict):
            data['templates'] = {}
        return data

    def _save_style_templates(self, data):
        """Persist classification style templates to JSON."""
        path = self._style_templates_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            json.dump(data, f, indent=4)
        return path

    def _normalise_category_style(self, style, fallback_index=0):
        """Return a complete, JSON-safe category style dictionary."""
        fallback = self._default_category_style(fallback_index)
        if not isinstance(style, dict):
            style = {}
        result = dict(fallback)
        result.update({k: v for k, v in style.items() if k in result})
        try:
            result['markersize'] = float(result.get('markersize', fallback['markersize']))
        except Exception:
            result['markersize'] = float(fallback['markersize'])
        try:
            result['linewidth'] = float(result.get('linewidth', fallback['linewidth']))
        except Exception:
            result['linewidth'] = float(fallback['linewidth'])
        try:
            result['alpha'] = float(result.get('alpha', fallback['alpha']))
        except Exception:
            result['alpha'] = float(fallback['alpha'])
        result['alpha'] = max(0.05, min(1.0, result['alpha']))
        for key in ('color', 'linecolor', 'arrowcolor', 'marker'):
            result[key] = str(result.get(key, fallback[key]))
        if result.get('fill') not in ('full', 'hollow'):
            result['fill'] = fallback['fill']
        return result

    def _set_artist_visible(self, artist, visible):
        """Set visibility on Matplotlib artists and ContourSet collections."""
        if artist is None:
            return
        if isinstance(artist, dict):
            return self._set_artist_visible(artist.get('artist'), visible)
        if isinstance(artist, (list, tuple)):
            for item in artist:
                self._set_artist_visible(item, visible)
            return
        if hasattr(artist, 'collections'):
            for item in artist.collections:
                item.set_visible(visible)
            return
        if hasattr(artist, 'set_visible'):
            artist.set_visible(visible)

    def _set_artist_alpha(self, artist, alpha):
        """Set alpha on a Matplotlib artist, or a dict/list of artists.

        Used to fade a category's raw points/lines when the interactive
        panel's "Show mean or Kamb maximum orientation" overlay is active, independently of
        the artist's stored category-style alpha (which _apply_category_style
        restores whenever a style is (re)applied).
        """
        if artist is None:
            return
        if isinstance(artist, dict):
            return self._set_artist_alpha(artist.get('artist'), alpha)
        if isinstance(artist, (list, tuple)):
            for item in artist:
                self._set_artist_alpha(item, alpha)
            return
        if hasattr(artist, 'set_alpha'):
            artist.set_alpha(alpha)

    def _bump_artist_zorder(self, artist, base_cache, bonus):
        """Offset a Matplotlib artist's (or dict/list/tuple's) z-order by
        `bonus`, relative to its own original z-order at construction time.

        The original z-order is cached in `base_cache` (keyed by id(artist))
        the first time each leaf artist is seen, so repeated calls with a
        different `bonus` always offset from the same baseline instead of
        compounding. This lets the interactive panel's category ordering
        control reorder which category's points/lines/mean-marker/label
        draw on top of another's, while preserving each layer's own
        raw/stats stacking tier (they get different base z-orders to begin
        with).
        """
        if artist is None:
            return
        if isinstance(artist, dict):
            return self._bump_artist_zorder(artist.get('artist'), base_cache, bonus)
        if isinstance(artist, (list, tuple)):
            for item in artist:
                self._bump_artist_zorder(item, base_cache, bonus)
            return
        if not hasattr(artist, 'set_zorder'):
            return
        key = id(artist)
        if key not in base_cache:
            base_cache[key] = artist.get_zorder() if hasattr(artist, 'get_zorder') else 1.0
        artist.set_zorder(base_cache[key] + bonus)

    def _open_category_panel(self, fig, artist_registry, category_counts=None,
                             category_styles=None, contour_update_callback=None,
                             girdle_update_callback=None,
                             title='Stereonet Categories', style_template_key=None,
                             export_legend_artists=None, show_visibility_controls=True,
                             stats_registry=None, label_registry=None, circle_registry=None):
        """Embed category visibility controls in a right-hand Qt panel.

        The panel is attached to the Matplotlib figure window as a Qt dock,
        rather than drawn on top of the stereonet axes. This prevents overlap
        with the stereonet, provides scrollbars for long category lists, and
        keeps the controls responsive when the figure window is resized.
        A Matplotlib-only fallback is kept for non-Qt backends.

        `stats_registry`/`label_registry`/`circle_registry` are the hidden
        per-category {'mean': ..., 'contour_max': ...} marker, text-label and
        (for planes) great-circle artists built by
        _build_orientation_stat_artists(). When `stats_registry` is
        non-empty, a "Statistics" control is added that lets the user pick
        which statistic to reveal (and, optionally, its label) and fades the
        raw points, without needing to re-plot.
        """
        if not artist_registry:
            return

        category_counts = category_counts or {}
        category_styles = category_styles or {}
        export_legend_artists = export_legend_artists or {}
        stats_registry = stats_registry or {}
        label_registry = label_registry or {}
        circle_registry = circle_registry or {}
        categories = sorted(artist_registry.keys(), key=lambda x: str(x))
        visible_state = {category: True for category in categories}
        category_order = list(categories)
        zorder_base_cache = {}
        order_spin_by_category = {}
        _stats_display_ref = [lambda: None]

        def _apply_category_order():
            for rank, category in enumerate(category_order):
                bonus = rank * 0.001
                for entry in artist_registry.get(category, []):
                    self._bump_artist_zorder(entry, zorder_base_cache, bonus)
                for entry in (stats_registry.get(category) or {}).values():
                    self._bump_artist_zorder(entry, zorder_base_cache, bonus)
                for label_entry in (label_registry.get(category) or {}).values():
                    self._bump_artist_zorder(label_entry, zorder_base_cache, bonus)
                for circle_entry in (circle_registry.get(category) or {}).values():
                    self._bump_artist_zorder(circle_entry, zorder_base_cache, bonus)
                for entry in export_legend_artists.get(category, []):
                    self._bump_artist_zorder(entry, zorder_base_cache, bonus)
            fig.canvas.draw_idle()

        def _reorder_category(category, new_rank_1indexed):
            new_index = max(0, min(len(category_order) - 1, new_rank_1indexed - 1))
            old_index = category_order.index(category)
            if new_index == old_index:
                return
            category_order.pop(old_index)
            category_order.insert(new_index, category)
            for cat in categories:
                spin = order_spin_by_category.get(cat)
                if spin is None:
                    continue
                spin.blockSignals(True)
                spin.setValue(category_order.index(cat) + 1)
                spin.blockSignals(False)
            _apply_category_order()

        def _mpl_colour_to_hex(colour, fallback='#000000'):
            try:
                return to_hex(to_rgba(colour or fallback))
            except Exception:
                return fallback

        def _apply_style_to_matplotlib_artist(artist, role, style):
            if artist is None:
                return
            if isinstance(artist, dict):
                return _apply_style_to_matplotlib_artist(
                    artist.get('artist'), artist.get('role', role), style)
            if isinstance(artist, (list, tuple)):
                for item in artist:
                    _apply_style_to_matplotlib_artist(item, role, style)
                return

            alpha = float(style.get('alpha', 1.0))
            marker_colour = style.get('color', '#000000')
            line_colour = style.get('linecolor', marker_colour)
            arrow_colour = style.get('arrowcolor', marker_colour)
            line_width = float(style.get('linewidth', 1.0))

            if role == 'legend_label':
                # Legend text should remain standard black text.  Category
                # style changes must update only the marker symbol, not the
                # label font colour.  Visibility is still handled separately
                # through _set_artist_visible.
                return

            if role == 'marker':
                hollow = style.get('fill') == 'hollow'
                marker_face = 'white' if hollow else marker_colour
                if hasattr(artist, 'set_marker'):
                    artist.set_marker(style.get('marker', 'o'))
                if hasattr(artist, 'set_markersize'):
                    artist.set_markersize(float(style.get('markersize', 5)))
                if hasattr(artist, 'set_markerfacecolor'):
                    artist.set_markerfacecolor(marker_face)
                if hasattr(artist, 'set_markeredgecolor'):
                    artist.set_markeredgecolor(marker_colour)
                if hasattr(artist, 'set_color'):
                    artist.set_color(marker_colour)
                if hasattr(artist, 'set_alpha'):
                    artist.set_alpha(alpha)
            elif role == 'arrow':
                if hasattr(artist, 'set_color'):
                    artist.set_color(arrow_colour)
                if hasattr(artist, 'set_edgecolor'):
                    artist.set_edgecolor(arrow_colour)
                if hasattr(artist, 'set_facecolor'):
                    artist.set_facecolor(arrow_colour)
                if hasattr(artist, 'set_linewidth'):
                    artist.set_linewidth(line_width)
                if hasattr(artist, 'set_alpha'):
                    artist.set_alpha(alpha)
            else:
                if hasattr(artist, 'set_color'):
                    artist.set_color(line_colour)
                if hasattr(artist, 'set_linewidth'):
                    artist.set_linewidth(line_width)
                if hasattr(artist, 'set_alpha'):
                    artist.set_alpha(alpha)

        def _apply_category_style(category):
            style = category_styles.get(category, self._default_category_style(0))
            for entry in artist_registry.get(category, []):
                if isinstance(entry, dict):
                    _apply_style_to_matplotlib_artist(entry.get('artist'), entry.get('role', 'marker'), style)
                else:
                    _apply_style_to_matplotlib_artist(entry, 'marker', style)
            for legend_artist in export_legend_artists.get(category, []):
                _apply_style_to_matplotlib_artist(legend_artist, 'marker', style)
            # Keep the mean/Kamb-maximum overlay's colour/marker/size (and,
            # for planes, its great circle's colour/width) in sync with the
            # individual points' style.
            _restyle_mean_artist(stats_registry.get(category), style)
            _restyle_mean_circle(circle_registry.get(category), style)
            # Re-applying a style resets the raw artist's alpha to the style's
            # own value, which would undo any "Show mean or Kamb maximum
            # orientation" fade - and the label's colour defaults to the
            # symbol colour - so reapply both on top, if the Statistics
            # toggle is on.
            _stats_display_ref[0]()
            fig.canvas.draw_idle()

        def _apply_all_category_styles():
            for category in categories:
                _apply_category_style(category)
                try:
                    style = category_styles.get(category, self._default_category_style(0))
                    if category in legend_symbol_by_category:
                        legend_symbol_by_category[category].set_symbol_style(
                            style.get('marker', 'o'), style.get('color', '#000000'),
                            hollow=style.get('fill') == 'hollow')
                        legend_symbol_by_category[category].setToolTip(
                            f"{style.get('marker', 'o')} / {style.get('color', '#000000')}")
                except NameError:
                    # The legend widgets are defined later in the Qt branch.
                    pass
            _refresh_category_visibility()

        def _refresh_category_visibility():
            for category, visible in visible_state.items():
                for artist in artist_registry.get(category, []):
                    self._set_artist_visible(artist, visible)
                for legend_artist in export_legend_artists.get(category, []):
                    self._set_artist_visible(legend_artist, visible)
            # The Statistics overlay's own visibility/fade also depends on
            # which categories are shown/hidden here.
            _stats_display_ref[0]()
            if contour_update_callback is not None:
                contour_update_callback(visible_state)
            if girdle_update_callback is not None:
                girdle_update_callback(visible_state)
            fig.canvas.draw_idle()

        manager = getattr(fig.canvas, 'manager', None)
        window = getattr(manager, 'window', None)

        # Preferred path: native Qt dock on the right-hand side of the
        # Matplotlib window. This avoids any overlay on the stereonet itself.
        if window is not None and hasattr(window, 'addDockWidget'):
            dock = QDockWidget(title, window)
            dock.setObjectName('StereonetCategoryDock')
            dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)

            panel = QWidget(dock)
            panel_layout = QVBoxLayout(panel)
            panel_layout.setContentsMargins(8, 8, 8, 8)
            panel_layout.setSpacing(6)

            if stats_registry:
                stats_group = QGroupBox('Statistics', panel)
                stats_layout = QVBoxLayout(stats_group)
                stats_caption = QLabel(
                    "Show each category's mean or Kamb maximum orientation, "
                    "with individual samples faded behind it.", stats_group)
                stats_caption.setWordWrap(True)
                stats_layout.addWidget(stats_caption)

                stats_checkbox = QCheckBox('Show mean or Kamb maximum orientation', stats_group)
                stats_layout.addWidget(stats_checkbox)

                mode_row = QHBoxLayout()
                mode_label = QLabel('Statistic:', stats_group)
                stat_mode_combo = QComboBox(stats_group)
                stat_mode_combo.addItem('Mean vector', 'mean')
                stat_mode_combo.addItem('Contour maximum (Kamb)', 'contour_max')
                stat_mode_combo.setEnabled(False)
                stat_mode_combo.setToolTip(
                    "'Mean vector' is the spherical (Fisher) mean. 'Contour "
                    "maximum' is the peak of a per-category Kamb density grid, "
                    "which can differ from the mean for multi-modal scatter.")
                mode_row.addWidget(mode_label)
                mode_row.addWidget(stat_mode_combo, 1)
                stats_layout.addLayout(mode_row)

                opacity_row = QHBoxLayout()
                opacity_label = QLabel('Background opacity:', stats_group)
                opacity_spin = QDoubleSpinBox(stats_group)
                opacity_spin.setRange(0.05, 1.0)
                opacity_spin.setDecimals(2)
                opacity_spin.setSingleStep(0.05)
                opacity_spin.setValue(0.30)
                opacity_spin.setEnabled(False)
                opacity_spin.setToolTip(
                    "Opacity of the individual points shown behind the "
                    "mean or Kamb maximum orientation overlay.")
                opacity_row.addWidget(opacity_label)
                opacity_row.addWidget(opacity_spin)
                stats_layout.addLayout(opacity_row)

                label_checkbox = QCheckBox('Show mean or Kamb maximum orientation as label', stats_group)
                label_checkbox.setEnabled(False)
                label_checkbox.setToolTip(
                    "Draw the mean or Kamb maximum orientation value as a text "
                    "label on the plot, instead of only showing the marker.")
                stats_layout.addWidget(label_checkbox)

                font_row = QHBoxLayout()
                font_label = QLabel('Label font:', stats_group)
                font_combo = QFontComboBox(stats_group)
                font_combo.setCurrentFont(QFont('Arial'))
                font_combo.setEnabled(False)
                font_row.addWidget(font_label)
                font_row.addWidget(font_combo, 1)
                stats_layout.addLayout(font_row)

                size_row = QHBoxLayout()
                size_label = QLabel('Label size:', stats_group)
                size_spin = QSpinBox(stats_group)
                size_spin.setRange(6, 24)
                size_spin.setValue(8)
                size_spin.setEnabled(False)
                size_row.addWidget(size_label)
                size_row.addWidget(size_spin)
                stats_layout.addLayout(size_row)

                panel_layout.addWidget(stats_group, 0)

                def _apply_stats_display():
                    enabled = stats_checkbox.isChecked()
                    show_label = label_checkbox.isChecked()
                    opacity = opacity_spin.value()
                    font_family = font_combo.currentFont().family()
                    font_size = size_spin.value()
                    mode = stat_mode_combo.currentData()
                    label_checkbox.setEnabled(enabled)
                    opacity_spin.setEnabled(enabled)
                    stat_mode_combo.setEnabled(enabled)
                    font_combo.setEnabled(enabled and show_label)
                    size_spin.setEnabled(enabled and show_label)
                    for category in categories:
                        cat_visible = visible_state.get(category, True)
                        style = category_styles.get(category, self._default_category_style(0))
                        marker_by_mode = stats_registry.get(category) or {}
                        for m, marker_artist in marker_by_mode.items():
                            self._set_artist_visible(marker_artist, enabled and cat_visible and m == mode)
                        circle_by_mode = circle_registry.get(category) or {}
                        for m, circle_artist in circle_by_mode.items():
                            self._set_artist_visible(circle_artist, enabled and cat_visible and m == mode)
                        label_by_mode = label_registry.get(category) or {}
                        for m, label_artist in label_by_mode.items():
                            label_artist.set_color(style.get('color', '#000000'))
                            label_artist.set_fontfamily(font_family)
                            label_artist.set_fontsize(font_size)
                            self._set_artist_visible(
                                label_artist, enabled and show_label and cat_visible and m == mode)
                        base_alpha = float(style.get('alpha', 1.0))
                        effective_alpha = opacity if enabled else base_alpha
                        for entry in artist_registry.get(category, []):
                            self._set_artist_alpha(entry, effective_alpha)
                    fig.canvas.draw_idle()

                _stats_display_ref[0] = _apply_stats_display
                stats_checkbox.stateChanged.connect(lambda _state: _apply_stats_display())
                stat_mode_combo.currentIndexChanged.connect(lambda _idx: _apply_stats_display())
                label_checkbox.stateChanged.connect(lambda _state: _apply_stats_display())
                opacity_spin.valueChanged.connect(lambda _val: _apply_stats_display())
                font_combo.currentFontChanged.connect(lambda _font: _apply_stats_display())
                size_spin.valueChanged.connect(lambda _val: _apply_stats_display())

            style_mgmt_group = QGroupBox('Style Management', panel)
            style_mgmt_layout = QVBoxLayout(style_mgmt_group)
            style_mgmt_layout.setContentsMargins(8, 8, 8, 8)
            style_mgmt_layout.setSpacing(6)
            style_mgmt_caption = QLabel(
                'Save, load, reset or delete reusable classification style templates.',
                style_mgmt_group)
            style_mgmt_caption.setWordWrap(True)
            style_mgmt_layout.addWidget(style_mgmt_caption)

            template_row = QHBoxLayout()
            save_template_btn = QPushButton('Save')
            load_template_btn = QPushButton('Load')
            reset_styles_btn = QPushButton('Reset')
            delete_template_btn = QPushButton('Delete')
            save_template_btn.setToolTip('Save the current category styles as a reusable template.')
            load_template_btn.setToolTip('Load a previously saved category style template.')
            reset_styles_btn.setToolTip('Reset category styles to the default palette and markers.')
            delete_template_btn.setToolTip('Delete an existing saved classification style template.')
            for btn in (save_template_btn, load_template_btn, reset_styles_btn, delete_template_btn):
                template_row.addWidget(btn)
            style_mgmt_layout.addLayout(template_row)
            panel_layout.addWidget(style_mgmt_group, 0)

            class_group = QGroupBox('Classification' if show_visibility_controls else 'Plot Styling', panel)
            class_layout = QVBoxLayout(class_group)
            class_layout.setContentsMargins(8, 8, 8, 8)
            class_layout.setSpacing(6)

            caption_text = (
                'Toggle category visibility, or set its drawing order (higher number = '
                'drawn on top). Contours are recalculated from visible records.'
                if show_visibility_controls
                else 'Edit the default plotting style for the current stereonet, or set '
                     'its drawing order (higher number = drawn on top).'
            )
            caption = QLabel(caption_text)
            caption.setWordWrap(True)
            class_layout.addWidget(caption)

            button_row = QHBoxLayout()
            show_btn = QPushButton('All')
            hide_btn = QPushButton('None')
            invert_btn = QPushButton('Invert')
            for btn in (show_btn, hide_btn, invert_btn):
                button_row.addWidget(btn)
            if show_visibility_controls:
                class_layout.addLayout(button_row)

            checkbox_by_category = {}

            # The category controls are scrollable.  The legend is deliberately
            # kept outside this scroll area and anchored at the bottom of the
            # dock, so it remains visually distinct from the on/off controls.
            controls_widget = QWidget(panel)
            controls_layout = QVBoxLayout(controls_widget)
            controls_layout.setContentsMargins(0, 0, 0, 0)
            controls_layout.setSpacing(4)

            style_button_by_category = {}

            def _style_dialog(category):
                style = category_styles.setdefault(category, self._default_category_style(0)).copy()
                dlg = QDialog(panel)
                dlg.setWindowTitle(f'Category style: {category}')
                layout = QVBoxLayout(dlg)

                form = QFormLayout()
                marker_cb = QComboBox(dlg)
                marker_options = [
                    ('Circle', 'o'), ('Square', 's'), ('Triangle up', '^'),
                    ('Triangle down', 'v'), ('Diamond', 'D'), ('Plus', 'P'),
                    ('Cross', 'X'), ('Star', '*'), ('Triangle left', '<'),
                    ('Triangle right', '>')]
                for label, marker in marker_options:
                    marker_cb.addItem(label, marker)
                marker_index = marker_cb.findData(style.get('marker', 'o'))
                if marker_index >= 0:
                    marker_cb.setCurrentIndex(marker_index)

                marker_size = QDoubleSpinBox(dlg)
                marker_size.setRange(1.0, 25.0)
                marker_size.setDecimals(1)
                marker_size.setSingleStep(0.5)
                marker_size.setValue(float(style.get('markersize', 5)))

                line_width = QDoubleSpinBox(dlg)
                line_width.setRange(0.1, 10.0)
                line_width.setDecimals(1)
                line_width.setSingleStep(0.2)
                line_width.setValue(float(style.get('linewidth', 1.0)))

                alpha = QDoubleSpinBox(dlg)
                alpha.setRange(0.05, 1.0)
                alpha.setDecimals(2)
                alpha.setSingleStep(0.05)
                alpha.setValue(float(style.get('alpha', 1.0)))

                colour_values = {
                    'color': _mpl_colour_to_hex(style.get('color', '#000000')),
                    'linecolor': _mpl_colour_to_hex(style.get('linecolor', style.get('color', '#000000'))),
                    'arrowcolor': _mpl_colour_to_hex(style.get('arrowcolor', style.get('color', '#000000'))),
                }

                def _colour_button(key, label):
                    btn = QPushButton(label, dlg)
                    btn.setStyleSheet(f'background-color: {colour_values[key]};')
                    def _choose_colour():
                        colour = QColorDialog.getColor(QColor(colour_values[key]), dlg, label)
                        if colour.isValid():
                            colour_values[key] = colour.name()
                            btn.setStyleSheet(f'background-color: {colour_values[key]};')
                    btn.clicked.connect(_choose_colour)
                    return btn

                marker_colour_btn = _colour_button('color', 'Symbol colour')
                line_colour_btn = _colour_button('linecolor', 'Line colour')
                arrow_colour_btn = _colour_button('arrowcolor', 'Arrow colour')

                hollow_checkbox = QCheckBox('Hollow (white fill, coloured outline)', dlg)
                hollow_checkbox.setChecked(style.get('fill') == 'hollow')

                form.addRow('Symbol shape:', marker_cb)
                form.addRow('Symbol size:', marker_size)
                form.addRow('Symbol colour:', marker_colour_btn)
                form.addRow('Symbol fill:', hollow_checkbox)
                form.addRow('Line colour:', line_colour_btn)
                form.addRow('Line width:', line_width)
                form.addRow('Transparency:', alpha)
                form.addRow('Arrow colour:', arrow_colour_btn)
                layout.addLayout(form)

                buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, dlg)
                layout.addWidget(buttons)
                buttons.accepted.connect(dlg.accept)
                buttons.rejected.connect(dlg.reject)

                if dlg.exec() == QDialog.Accepted:
                    style.update({
                        'marker': marker_cb.currentData(),
                        'markersize': float(marker_size.value()),
                        'color': colour_values['color'],
                        'linecolor': colour_values['linecolor'],
                        'arrowcolor': colour_values['arrowcolor'],
                        'linewidth': float(line_width.value()),
                        'alpha': float(alpha.value()),
                        'fill': 'hollow' if hollow_checkbox.isChecked() else 'full',
                    })
                    category_styles[category] = style
                    _apply_category_style(category)
                    if category in legend_symbol_by_category:
                        legend_symbol_by_category[category].set_symbol_style(
                            style.get('marker', 'o'), style.get('color', '#000000'),
                            hollow=style.get('fill') == 'hollow')
                        legend_symbol_by_category[category].setToolTip(
                            f"{style.get('marker', 'o')} / {style.get('color', '#000000')}")
                    _refresh_category_visibility()

            for category in categories:
                style = category_styles.get(category, self._default_category_style(0))
                row_widget = QWidget(controls_widget)
                row_layout = QHBoxLayout(row_widget)
                row_layout.setContentsMargins(0, 0, 0, 0)
                row_layout.setSpacing(4)

                checkbox = QCheckBox(f'{category} (n={category_counts.get(category, 0)})', row_widget)
                checkbox.setChecked(True)
                checkbox.setVisible(show_visibility_controls)
                checkbox_by_category[category] = checkbox

                order_spin = QSpinBox(row_widget)
                order_spin.setRange(1, len(categories))
                order_spin.setValue(category_order.index(category) + 1)
                order_spin.setToolTip(
                    'Drawing (z-stack) order: higher numbers draw on top of lower ones.')
                order_spin.setFixedWidth(44)
                order_spin_by_category[category] = order_spin

                style_btn = QPushButton('Style…', row_widget)
                style_btn.setToolTip(f'Edit plotting style for {category}')
                style_button_by_category[category] = style_btn

                def _make_state_callback(cat):
                    def _on_state_changed(state):
                        visible_state[cat] = (state == Qt.Checked)
                        _refresh_category_visibility()
                    return _on_state_changed

                def _make_style_callback(cat):
                    return lambda: _style_dialog(cat)

                def _make_order_callback(cat):
                    return lambda value: _reorder_category(cat, value)

                checkbox.stateChanged.connect(_make_state_callback(category))
                style_btn.clicked.connect(_make_style_callback(category))
                order_spin.valueChanged.connect(_make_order_callback(category))
                if show_visibility_controls:
                    row_layout.addWidget(checkbox, 1)
                else:
                    row_layout.addWidget(QLabel('Default style', row_widget), 1)
                row_layout.addWidget(order_spin, 0)
                row_layout.addWidget(style_btn, 0)
                controls_layout.addWidget(row_widget)

            controls_layout.addStretch(1)

            scroll = QScrollArea(dock)
            scroll.setWidgetResizable(True)
            scroll.setWidget(controls_widget)
            scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            scroll.setMinimumHeight(120)
            class_layout.addWidget(scroll, 1)
            panel_layout.addWidget(class_group, 1)

            class _MarkerSymbolWidget(QWidget):
                """Vector marker preview for the classification legend.

                The symbol is painted directly with QPainter instead of being
                rendered to a raster pixmap. This keeps the legend crisp when
                the dock or plot window is resized and uses the same marker
                shape/colour assigned to the plotted category.
                """

                def __init__(self, marker, colour, parent=None, hollow=False):
                    super().__init__(parent)
                    self.marker = marker or 'o'
                    self.hollow = bool(hollow)
                    # QColor does not understand Matplotlib colour names such as
                    # ``tab:orange`` or ``tab:blue``. Convert every Matplotlib-
                    # compatible colour to a hex string first, then fall back to
                    # black only if the value is genuinely invalid.
                    try:
                        qcolour = QColor(to_hex(to_rgba(colour or '#000000')))
                    except Exception:
                        qcolour = QColor('#000000')
                    self.colour = qcolour
                    self.setMinimumSize(24, 24)
                    self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

                def sizeHint(self):
                    return QSize(24, 24)

                def set_symbol_style(self, marker, colour, hollow=False):
                    self.marker = marker or 'o'
                    self.hollow = bool(hollow)
                    try:
                        self.colour = QColor(to_hex(to_rgba(colour or '#000000')))
                    except Exception:
                        self.colour = QColor('#000000')
                    self.update()

                def paintEvent(self, event):
                    painter = QPainter(self)
                    painter.setRenderHint(QPainter.Antialiasing)

                    side = min(self.width(), self.height())
                    cx = self.width() / 2.0
                    cy = self.height() / 2.0
                    r = side * 0.30

                    fill_colour = QColor('white') if self.hollow else self.colour
                    pen = QPen(self.colour)
                    pen.setWidthF(max(1.2, side * 0.09))
                    painter.setPen(pen)
                    painter.setBrush(QBrush(fill_colour))

                    m = self.marker
                    if m == 'o':
                        painter.drawEllipse(QPointF(cx, cy), r, r)
                    elif m == 's':
                        painter.drawRect(QRectF(cx - r, cy - r, 2 * r, 2 * r))
                    elif m == '^':
                        painter.drawPolygon(QPolygonF([
                            QPointF(cx, cy - r),
                            QPointF(cx - r, cy + r),
                            QPointF(cx + r, cy + r),
                        ]))
                    elif m == 'v':
                        painter.drawPolygon(QPolygonF([
                            QPointF(cx - r, cy - r),
                            QPointF(cx + r, cy - r),
                            QPointF(cx, cy + r),
                        ]))
                    elif m == '<':
                        painter.drawPolygon(QPolygonF([
                            QPointF(cx - r, cy),
                            QPointF(cx + r, cy - r),
                            QPointF(cx + r, cy + r),
                        ]))
                    elif m == '>':
                        painter.drawPolygon(QPolygonF([
                            QPointF(cx + r, cy),
                            QPointF(cx - r, cy - r),
                            QPointF(cx - r, cy + r),
                        ]))
                    elif m == 'D':
                        painter.drawPolygon(QPolygonF([
                            QPointF(cx, cy - r),
                            QPointF(cx + r, cy),
                            QPointF(cx, cy + r),
                            QPointF(cx - r, cy),
                        ]))
                    elif m in ('P', '+'):
                        # Plus/cross-style markers are stroked with the category
                        # colour.  ``P`` also gets a small filled centre to match
                        # Matplotlib's filled-plus marker more closely.
                        painter.setBrush(Qt.NoBrush)
                        painter.drawLine(QPointF(cx - r, cy), QPointF(cx + r, cy))
                        painter.drawLine(QPointF(cx, cy - r), QPointF(cx, cy + r))
                        if m == 'P':
                            painter.setBrush(QBrush(fill_colour))
                            painter.drawRect(QRectF(cx - r * 0.38, cy - r * 0.38,
                                                    r * 0.76, r * 0.76))
                    elif m in ('X', 'x'):
                        painter.setBrush(Qt.NoBrush)
                        painter.drawLine(QPointF(cx - r, cy - r), QPointF(cx + r, cy + r))
                        painter.drawLine(QPointF(cx - r, cy + r), QPointF(cx + r, cy - r))
                    elif m == '*':
                        painter.setBrush(Qt.NoBrush)
                        painter.drawLine(QPointF(cx - r, cy), QPointF(cx + r, cy))
                        painter.drawLine(QPointF(cx, cy - r), QPointF(cx, cy + r))
                        painter.drawLine(QPointF(cx - r * 0.72, cy - r * 0.72),
                                         QPointF(cx + r * 0.72, cy + r * 0.72))
                        painter.drawLine(QPointF(cx - r * 0.72, cy + r * 0.72),
                                         QPointF(cx + r * 0.72, cy - r * 0.72))
                    else:
                        painter.drawEllipse(QPointF(cx, cy), r, r)

                    painter.end()

            legend_group = QGroupBox('Legend', panel)
            legend_symbol_by_category = {}
            legend_layout = QVBoxLayout(legend_group)
            legend_layout.setContentsMargins(8, 8, 8, 8)
            legend_layout.setSpacing(4)

            for category in categories:
                style = category_styles.get(category, self._default_category_style(0))
                row_widget = QWidget(legend_group)
                row_layout = QHBoxLayout(row_widget)
                row_layout.setContentsMargins(0, 0, 0, 0)
                row_layout.setSpacing(6)

                marker = style.get('marker', 'o')
                colour = style.get('color', '#000000')
                marker_label = _MarkerSymbolWidget(
                    marker, colour, row_widget, hollow=style.get('fill') == 'hollow')
                marker_label.setToolTip(f'{marker} / {colour}')

                category_label = QLabel(f'{category} (n={category_counts.get(category, 0)})', row_widget)
                category_label.setWordWrap(True)

                row_layout.addWidget(marker_label)
                row_layout.addWidget(category_label, 1)
                legend_layout.addWidget(row_widget)
                legend_symbol_by_category[category] = marker_label

            legend_scroll = QScrollArea(dock)
            legend_scroll.setWidgetResizable(True)
            legend_scroll.setWidget(legend_group)
            legend_scroll.setMaximumHeight(180)
            panel_layout.addWidget(legend_scroll, 0)

            dock.setWidget(panel)
            dock.setMinimumWidth(260)
            dock.setFeatures(QDockWidget.DockWidgetMovable |
                             QDockWidget.DockWidgetFloatable)
            window.addDockWidget(Qt.RightDockWidgetArea, dock)

            def _set_checkbox_state(category, state):
                checkbox = checkbox_by_category[category]
                if checkbox.isChecked() != state:
                    checkbox.blockSignals(True)
                    checkbox.setChecked(state)
                    checkbox.blockSignals(False)

            def _show_all():
                for category in categories:
                    visible_state[category] = True
                    _set_checkbox_state(category, True)
                _refresh_category_visibility()

            def _hide_all():
                for category in categories:
                    visible_state[category] = False
                    _set_checkbox_state(category, False)
                _refresh_category_visibility()

            def _invert():
                for category in categories:
                    visible_state[category] = not visible_state[category]
                    _set_checkbox_state(category, visible_state[category])
                _refresh_category_visibility()

            def _template_key():
                key = style_template_key or 'default'
                return str(key) if str(key).strip() else 'default'

            def _template_name(default_name=None):
                default_name = default_name or _template_key()
                name, ok = QInputDialog.getText(
                    panel, 'Classification style template',
                    'Template name:', text=str(default_name))
                if not ok:
                    return None
                name = str(name).strip()
                return name or None

            def _sync_legend_symbols():
                for category in categories:
                    style = category_styles.get(category, self._default_category_style(0))
                    if category in legend_symbol_by_category:
                        legend_symbol_by_category[category].set_symbol_style(
                            style.get('marker', 'o'), style.get('color', '#000000'),
                            hollow=style.get('fill') == 'hollow')
                        legend_symbol_by_category[category].setToolTip(
                            f"{style.get('marker', 'o')} / {style.get('color', '#000000')}")

            def _apply_all_category_styles():
                for category in categories:
                    _apply_category_style(category)
                _sync_legend_symbols()
                _refresh_category_visibility()

            def _save_template():
                name = _template_name(_template_key())
                if not name:
                    return
                data = self._load_style_templates()
                templates = data.setdefault('templates', {})
                templates[name] = {
                    'classificationField': _template_key(),
                    'styles': {
                        str(category): self._normalise_category_style(
                            category_styles.get(category, self._default_category_style(i)), i)
                        for i, category in enumerate(categories)
                    }
                }
                try:
                    path = self._save_style_templates(data)
                except Exception as exc:
                    QMessageBox.critical(panel, 'Save styles',
                                         f'Could not save style template:\n{exc}')
                    return
                QMessageBox.information(panel, 'Save styles',
                                        f'Style template "{name}" saved to:\n{path}')

            def _load_template():
                data = self._load_style_templates()
                templates = data.get('templates', {}) if isinstance(data, dict) else {}
                if not templates:
                    QMessageBox.information(panel, 'Load styles',
                                            'No saved classification style templates were found.')
                    return
                names = sorted(templates.keys(), key=lambda x: str(x).lower())
                preferred = _template_key()
                current_index = names.index(preferred) if preferred in names else 0
                name, ok = QInputDialog.getItem(
                    panel, 'Load styles', 'Select a style template:',
                    names, current_index, False)
                if not ok or not name:
                    return
                template = templates.get(str(name), {})
                styles = template.get('styles', {}) if isinstance(template, dict) else {}
                if not isinstance(styles, dict) or not styles:
                    QMessageBox.warning(panel, 'Load styles',
                                        'The selected style template does not contain any category styles.')
                    return
                for i, category in enumerate(categories):
                    if str(category) in styles:
                        category_styles[category] = self._normalise_category_style(styles[str(category)], i)
                _apply_all_category_styles()

            def _reset_styles():
                reply = QMessageBox.question(
                    panel, 'Reset styles',
                    'Reset all category styles to the default palette?',
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                if reply != QMessageBox.Yes:
                    return
                for i, category in enumerate(categories):
                    category_styles[category] = self._default_category_style(i)
                _apply_all_category_styles()

            def _delete_template():
                data = self._load_style_templates()
                templates = data.get('templates', {}) if isinstance(data, dict) else {}
                if not templates:
                    QMessageBox.information(panel, 'Delete styles',
                                            'No saved classification style templates were found.')
                    return
                names = sorted(templates.keys(), key=lambda x: str(x).lower())
                preferred = _template_key()
                current_index = names.index(preferred) if preferred in names else 0
                name, ok = QInputDialog.getItem(
                    panel, 'Delete styles', 'Select a style template to delete:',
                    names, current_index, False)
                if not ok or not name:
                    return
                reply = QMessageBox.question(
                    panel, 'Delete styles',
                    f'Delete style template "{name}" permanently?',
                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                if reply != QMessageBox.Yes:
                    return
                templates.pop(str(name), None)
                data['templates'] = templates
                try:
                    path = self._save_style_templates(data)
                except Exception as exc:
                    QMessageBox.critical(panel, 'Delete styles',
                                         f'Could not delete style template:\n{exc}')
                    return
                QMessageBox.information(panel, 'Delete styles',
                                        f'Style template "{name}" deleted from:\n{path}')

            show_btn.clicked.connect(_show_all)
            hide_btn.clicked.connect(_hide_all)
            invert_btn.clicked.connect(_invert)
            save_template_btn.clicked.connect(_save_template)
            load_template_btn.clicked.connect(_load_template)
            reset_styles_btn.clicked.connect(_reset_styles)
            delete_template_btn.clicked.connect(_delete_template)

            self._category_controls = {
                'dock': dock,
                'panel': panel,
                'scroll': scroll,
                'legend': legend_group,
                'legend_scroll': legend_scroll,
                'checkboxes': checkbox_by_category,
                'style_buttons': style_button_by_category,
                'styles': category_styles,
                'buttons': (show_btn, hide_btn, invert_btn),
                'style_management_group': style_mgmt_group,
                'classification_group': class_group,
                'style_template_buttons': (save_template_btn, load_template_btn, reset_styles_btn, delete_template_btn),
                'visible_state': visible_state,
            }

            if contour_update_callback is not None:
                contour_update_callback(visible_state)
            if girdle_update_callback is not None:
                girdle_update_callback(visible_state)
            return

        # Fallback for non-Qt Matplotlib backends. Reserve space outside the
        # stereonet axes and draw controls in that reserved figure margin.
        labels = [f'{category} (n={category_counts.get(category, 0)})'
                  for category in categories]
        label_to_category = dict(zip(labels, categories))
        try:
            fig.subplots_adjust(right=0.68)
        except Exception:
            pass

        check_ax = fig.add_axes([0.72, 0.42, 0.25, 0.45])
        check_ax.set_title('Categories', fontsize=9)
        checks = CheckButtons(check_ax, labels, [True] * len(labels))

        def _set_button_state(index, state):
            if checks.get_status()[index] != state:
                checks.set_active(index)

        def _on_clicked(label):
            category = label_to_category.get(label)
            if category is None:
                return
            index = categories.index(category)
            visible_state[category] = checks.get_status()[index]
            _refresh_category_visibility()

        checks.on_clicked(_on_clicked)
        show_ax = fig.add_axes([0.72, 0.32, 0.075, 0.05])
        hide_ax = fig.add_axes([0.81, 0.32, 0.075, 0.05])
        invert_ax = fig.add_axes([0.90, 0.32, 0.075, 0.05])
        show_btn = Button(show_ax, 'All')
        hide_btn = Button(hide_ax, 'None')
        invert_btn = Button(invert_ax, 'Invert')

        def _show_all(event=None):
            for i, category in enumerate(categories):
                visible_state[category] = True
                _set_button_state(i, True)
            _refresh_category_visibility()

        def _hide_all(event=None):
            for i, category in enumerate(categories):
                visible_state[category] = False
                _set_button_state(i, False)
            _refresh_category_visibility()

        def _invert(event=None):
            for i, category in enumerate(categories):
                visible_state[category] = not visible_state[category]
                _set_button_state(i, visible_state[category])
            _refresh_category_visibility()

        show_btn.on_clicked(_show_all)
        hide_btn.on_clicked(_hide_all)
        invert_btn.on_clicked(_invert)
        self._category_controls = {
            'checkbuttons': checks,
            'buttons': (show_btn, hide_btn, invert_btn),
            'visible_state': visible_state,
            'axes': (check_ax, show_ax, hide_ax, invert_ax),
        }
        if contour_update_callback is not None:
            contour_update_callback(visible_state)
        if girdle_update_callback is not None:
            girdle_update_callback(visible_state)

    def rose_diagram(self, strikes, title):
        """Plot an azimuth rose diagram with the same visual footprint as the stereonet.

        The polar axes are deliberately smaller than the figure canvas so the
        layer title and azimuth tick labels do not overlap.
        """
        bin_edges = np.arange(-5, 366, 10)
        number_of_strikes, bin_edges = np.histogram(strikes, bin_edges)
        number_of_strikes[0] += number_of_strikes[-1]
        half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
        two_halves = np.concatenate([half, half])

        fig = plt.figure(figsize=(5.4, 5.4))
        ax = fig.add_subplot(111, projection='polar')
        # Keep the rose itself smaller than the canvas, comparable to the
        # stereonet circle and clear of the title/azimuth labels.
        ax.set_position([0.18, 0.12, 0.64, 0.64])

        ax.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves,
               width=np.deg2rad(10), color='.8', edgecolor='k')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10), fontsize=7)
        if two_halves.max() > 0:
            step = max(1, int(np.ceil(two_halves.max() / 5.0)))
            ax.set_rgrids(np.arange(step, two_halves.max() + 1, step), angle=75, fontsize=7)
            for _rose_label in ax.get_xticklabels() + ax.get_yticklabels():
                _rose_label.set_fontsize(7)
        ax.set_title(title, y=1.22, fontsize=10, pad=10)
        plt.show()
        
    def contourPlot(self):
        names = self._structural_field_names()
        snames = names['strike']
        ddnames = names['dipdir']
        dnames = names['dip']
        anames = names['azimuth']
        pnames = names['plunge']
        srefnames = names['strike_ref']
        drefnames = names['dip_ref']
        ddrefnames = names['dipdir_ref']
        knames = names['kinematics']
        prhrnames = names['pitch_rhr']
        
        plane_strikes = list()
        plane_dips = list()
        plane_feature_ids = list()
        plane_labels = list()
        plane_categories = list()
        linear_plunges = list()
        linear_bearings = list()
        linear_feature_ids = list()
        linear_labels = list()
        linear_categories = list()
        strikesref = list()
        dipsref = list()
        ref_categories = list()
        plunges = list()
        kinematics = list()
        rakes_strikes = list()
        rakes_dips = list()
        roseAzimuth = list()
        rhr = list()
        azs = list()
        kinematic_arrow_records = list()


        project = QgsProject.instance()
        proj_file_path=project.fileName()
        head_tail = os.path.split(proj_file_path)
        WAXI_project_path = os.path.abspath(QgsProject.instance().fileName())
        stereoConfigPath = os.path.join(os.path.dirname(WAXI_project_path), "99_COMMAND_FILES_PLUGIN/stereonet.json")

        #stereoConfigPath = head_tail[0]+"/0. FIELD DATA/0. CURRENT MISSION/0. STOPS-SAMPLING-PHOTOGRAPHS-COMMENTS/stereonet.json"
        
        stereoConfig = {'showGtCircles': False, 'showContours': True,
                        'showKinematics': False, 'linPlanes': True, 'roseDiagram': False,
                        'fitGirdle': False, 'dataType': 'Planes Only',
                        'kinematicsField': None,
                        'kinematicsAnchor': 'Plane pole',
                        'classificationEnabled': False,
                        'classificationField': None,
                        'filterEnabled': False,
                        'filterExpression': ''}

        if os.path.exists(stereoConfigPath):
            with open(stereoConfigPath, "r") as json_file:
                stereoConfig = json.load(json_file)
        else:
            qs = StereonetSettingsDialog.load_qsettings()
            if qs is not None:
                stereoConfig = qs
        
        self.iface.layerTreeView().selectedLayers()

        layers = self.iface.layerTreeView().selectedLayers()
        active_layer = self.iface.activeLayer()
        if (not layers and active_layer is not None and
                active_layer.type() == QgsMapLayer.VectorLayer):
            layers = [active_layer]
        elif active_layer is not None and active_layer.type() == QgsMapLayer.VectorLayer:
            # If the layer-tree selection is stale but the active layer carries
            # the selected structural features, prefer the active layer.  This
            # avoids the misleading "No data selected" warning after switching
            # between layers without reopening the settings dialog.
            try:
                selected_counts = [len(layer.selectedFeatures()) for layer in layers
                                   if layer.type() == QgsMapLayer.VectorLayer]
                if selected_counts and sum(selected_counts) == 0 and len(active_layer.selectedFeatures()) > 0:
                    layers = [active_layer]
            except Exception:
                pass

        current_plot_layer_signature = '|'.join(sorted(
            str(layer.id()) for layer in layers
            if layer.type() == QgsMapLayer.VectorLayer))
        saved_plot_layer_signature = str(
            stereoConfig.get('plotLayerSignature') or
            stereoConfig.get('classificationLayerSignature', '') or '')
        if (stereoConfig.get('roseDiagram', False) and current_plot_layer_signature and
                current_plot_layer_signature != saved_plot_layer_signature):
            stereoConfig['roseDiagram'] = False

        selected_kinematics_field = stereoConfig.get('kinematicsField')
        if stereoConfig.get('showKinematics', False):
            # Saved kinematics settings can become stale when switching layers.
            # If the saved field is absent, automatically fall back to the first
            # valid kinematics candidate in the currently selected/active layer.
            field_available = False
            if selected_kinematics_field:
                for _layer in layers:
                    if (_layer.type() == QgsMapLayer.VectorLayer and
                            _layer.fields().lookupField(selected_kinematics_field) != -1):
                        field_available = True
                        break
            if not field_available:
                _candidates = self._candidate_kinematics_fields(layers)
                selected_kinematics_field = _candidates[0] if _candidates else None
                stereoConfig['kinematicsField'] = selected_kinematics_field

        plot_kinematics = bool(stereoConfig.get('showKinematics', False)
                               and selected_kinematics_field)

        # Resolve the effective plotting mode before extracting feature
        # orientations.  The extraction loop needs this value for manual
        # overrides such as plotting DipDir/Dip as Trend/Plunge in Lineations
        # Only mode.
        detected_data_type = self._detect_data_type_from_layers(layers)
        saved_data_type = stereoConfig.get('dataType', '') or ''
        if saved_data_type == 'Lineations with Planes':
            saved_data_type = 'Lineations with Bearing Planes'
        valid_data_types = {
            'Planes Only',
            'Lineations Only',
            'Lineations with Bearing Planes'
        }
        effective_data_type = saved_data_type if saved_data_type in valid_data_types else detected_data_type
        effective_data_type = effective_data_type or 'Planes Only'
        classification_enabled = bool(stereoConfig.get('classificationEnabled', False))
        classification_field = stereoConfig.get('classificationField') if classification_enabled else None

        if classification_enabled and classification_field:
            # Saved classification settings can become stale when changing
            # layer.  If the field is absent or only contains NULL/empty values
            # in the currently selected features, silently fall back to an
            # unclassified plot instead of discarding every feature.
            has_class_value = False
            for _layer in layers:
                if _layer.type() != QgsMapLayer.VectorLayer:
                    continue
                if _layer.fields().lookupField(classification_field) == -1:
                    continue
                for _feature in _layer.selectedFeatures():
                    _value = _attr(_feature[classification_field])
                    if _value is not None and str(_value).strip() != '':
                        has_class_value = True
                        break
                if has_class_value:
                    break
            if not has_class_value:
                classification_enabled = False
                classification_field = None
                stereoConfig['classificationEnabled'] = False
                stereoConfig['classificationField'] = None

        filter_expression = self._build_filter_expression(
            stereoConfig.get('filterExpression', '') if stereoConfig.get('filterEnabled', False) else '')
        if filter_expression is not None:
            try:
                referenced = set(filter_expression.referencedColumns())
            except Exception:
                referenced = set()
            if referenced:
                available = set()
                for _layer in layers:
                    if _layer.type() == QgsMapLayer.VectorLayer:
                        available.update([field.name() for field in _layer.fields()])
                missing = [field for field in referenced if field not in available]
                if missing:
                    self.iface.messageBar().pushMessage(
                        'Stereonet filter ignored',
                        'The saved filter references field(s) not present in the selected layer: ' +
                        ', '.join(missing),
                        level=Qgis.Warning, duration=8)
                    filter_expression = None

        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer:

                iter = layer.selectedFeatures()
                if filter_expression is not None:
                    iter = [feature for feature in iter
                            if self._feature_passes_filter(layer, feature, filter_expression)]
                class_field_exists = bool(classification_field and
                                          layer.fields().lookupField(classification_field) != -1)
                strikeExists, sname = self._field_exists(layer,snames)
                ddrExists, ddname = self._field_exists(layer,ddnames)
                dipExists, dname = self._field_exists(layer,dnames)
                azimuthExists, aname = self._field_exists(layer,anames)
                plungeExists, pname = self._field_exists(layer,pnames)
                srefExists, srefname = self._field_exists(layer,srefnames)
                drefExists, drefname = self._field_exists(layer,drefnames)
                ddrefExists, ddrefname = self._field_exists(layer,ddrefnames)
                if (plot_kinematics and selected_kinematics_field and
                        layer.fields().lookupField(selected_kinematics_field) != -1):
                    kinematicsExists, kname = True, selected_kinematics_field
                else:
                    kinematicsExists, kname = self._field_exists(layer,knames)
                prhrExists, prhrname = self._field_exists(layer,prhrnames )

                # If the user deliberately forces a planar layer to plot as
                # Lineations Only and there are no Trend/Azimuth + Plunge
                # fields, use DipDir as Trend and Dip as Plunge.  This is a
                # manual-override fallback only; native lineation fields still
                # take priority whenever they exist.
                use_plane_orientation_as_lineation = (
                    effective_data_type == 'Lineations Only' and
                    not (azimuthExists and plungeExists) and
                    ddrExists and dipExists
                )



                for feature in iter:
                    current_plane_strike = None
                    current_plane_dip = None
                    current_ref_strike = None
                    current_ref_dip = None
                    current_line_plunge = None
                    current_line_bearing = None
                    if classification_enabled and class_field_exists:
                        current_category = self._category_value(feature, classification_field)
                    else:
                        current_category = 'All'
                    if current_category is None:
                        continue

                    # Capture plane data (dip direction/dip or strike/dip)
                    if ddrExists and dipExists:
                        val_dd, val_d = _attr(feature[ddname]), _attr(feature[dname])
                        if val_dd is not None and val_d is not None:
                            current_plane_strike = (float(val_dd) - 90.0) % 360.0
                            current_plane_dip = float(val_d)
                            plane_strikes.append(current_plane_strike)
                            plane_dips.append(current_plane_dip)
                            plane_feature_ids.append((layer, feature.id()))
                            plane_labels.append(f"{int(val_d)}/{int(val_dd):03d}")
                            plane_categories.append(current_category)
                    elif strikeExists and dipExists:
                        val_s, val_d = _attr(feature[sname]), _attr(feature[dname])
                        if val_s is not None and val_d is not None:
                            current_plane_strike = float(val_s)
                            current_plane_dip = float(val_d)
                            plane_strikes.append(current_plane_strike)
                            plane_dips.append(current_plane_dip)
                            plane_feature_ids.append((layer, feature.id()))
                            plane_labels.append(f"{int(val_d)}/{int(val_s):03d}")
                            plane_categories.append(current_category)

                    # Capture linear data (azimuth/plunge) independently.
                    # Manual fallback: if Lineations Only is selected on a
                    # planar dataset with no Trend/Plunge fields, use
                    # DipDir as Trend and Dip as Plunge.
                    if azimuthExists and plungeExists:
                        val_a, val_p = _attr(feature[aname]), _attr(feature[pname])
                    elif use_plane_orientation_as_lineation:
                        val_a, val_p = _attr(feature[ddname]), _attr(feature[dname])
                    else:
                        val_a, val_p = None, None

                    if val_a is not None and val_p is not None:
                        current_line_plunge = float(val_p)
                        current_line_bearing = float(val_a)
                        linear_plunges.append(current_line_plunge)
                        linear_bearings.append(current_line_bearing)
                        linear_feature_ids.append((layer, feature.id()))
                        linear_labels.append(f"{int(float(val_p))}/{int(float(val_a)):03d}")
                        linear_categories.append(current_category)

                    if drefExists:
                        vd = _attr(feature[drefname])
                        if vd is not None and srefExists:
                            vs = _attr(feature[srefname])
                            if vs is not None:
                                current_ref_strike = float(vs)
                                current_ref_dip = float(vd)
                                strikesref.append(current_ref_strike)
                                dipsref.append(current_ref_dip)
                                ref_categories.append(current_category)
                        elif vd is not None and ddrefExists:
                            vdd = _attr(feature[ddrefname])
                            if vdd is not None:
                                current_ref_strike = (float(vdd) - 90.0) % 360.0
                                current_ref_dip = float(vd)
                                strikesref.append(current_ref_strike)
                                dipsref.append(current_ref_dip)
                                ref_categories.append(current_category)

                    if (plungeExists and drefExists and
                            kinematicsExists and azimuthExists):
                        vp = _attr(feature[pname])
                        vdr = _attr(feature[drefname])
                        vk = _attr(feature[kname])
                        if vp is not None and vdr is not None and vk is not None:
                            if srefExists:
                                rake_strike = _attr(feature[srefname])
                            elif ddrefExists and _attr(feature[ddrefname]) is not None:
                                rake_strike = (float(_attr(feature[ddrefname])) - 90) % 360
                            else:
                                rake_strike = None
                            if rake_strike is not None:
                                rakes_strikes.append(rake_strike)
                                rakes_dips.append(vdr)
                                kinematics.append(vk)
                                rhr.append(_attr(feature[prhrname]) if prhrExists else None)
                                azs.append(_attr(feature[aname]))

                    if plot_kinematics and kinematicsExists and current_line_bearing is not None and current_line_plunge is not None:
                        vk = _attr(feature[kname])
                        arrow_strike = current_ref_strike
                        arrow_dip = current_ref_dip
                        if arrow_strike is None or arrow_dip is None:
                            arrow_strike = current_plane_strike
                            arrow_dip = current_plane_dip

                        self._append_kinematic_arrow_record(
                            kinematic_arrow_records,
                            arrow_strike, arrow_dip,
                            current_line_plunge, current_line_bearing,
                            vk, current_category)

                    if stereoConfig.get('roseDiagram', False):
                        # The rose diagram should use the orientation actually
                        # relevant to the plotted data. For lineation-bearing
                        # layers this is the lineation bearing/trend; for
                        # plane-only layers this falls back to the plane strike
                        # derived from Strike or Dip Direction. The previous
                        # implementation only read Azimuth/Trend fields, so
                        # planar datasets produced an empty rose diagram.
                        if current_line_bearing is not None:
                            roseAzimuth.append(current_line_bearing % 360.0)
                        elif current_plane_strike is not None:
                            roseAzimuth.append(current_plane_strike % 360.0)
 


            else:
                continue

        ref_clean = [(sref, dref, cat) for sref, dref, cat in zip(strikesref, dipsref, ref_categories)
                     if sref is not None and dref is not None]
        if ref_clean:
            strikesref, dipsref, ref_categories = map(list, zip(*ref_clean))
        else:
            strikesref, dipsref, ref_categories = [], [], []

        show_planes = effective_data_type == 'Planes Only'
        show_linears = effective_data_type in (
            'Lineations Only',
            'Lineations with Bearing Planes')
        show_bearing_planes = effective_data_type == 'Lineations with Bearing Planes'
        # In combined lineation/reference-plane mode, the bearing planes are
        # part of the selected data type itself.  Do not rely on a stale saved
        # linPlanes value from QSettings/stereonet.json; otherwise the checkbox
        # can appear ON in the dialog while the first plot still suppresses the
        # reference-plane great circles until the user toggles and re-saves it.
        effective_lin_planes = (
            True if show_bearing_planes
            else stereoConfig.get('linPlanes', True)
        )
        has_planes = len(plane_strikes) > 0
        has_linears = len(linear_plunges) > 0
        has_bearing_planes = len(strikesref) > 0 and len(dipsref) > 0

        if len(roseAzimuth) != 0 and stereoConfig.get('roseDiagram', False):
            rose_layers = [l.name() for l in layers if l.type() == QgsMapLayer.VectorLayer]
            rose_title = (', '.join(rose_layers) if rose_layers else 'Selected data')
            self.rose_diagram(roseAzimuth, rose_title + " [# " + str(len(roseAzimuth)) + "]")
        elif ((show_planes and has_planes) or
              (show_linears and has_linears) or
              (show_bearing_planes and has_bearing_planes)):
            fig, ax = mplstereonet.subplots()
            # Reserve a narrow right-hand margin for the exportable legend.
            # The Qt category dock remains outside the saved figure.
            try:
                ax.set_position([0.24, 0.10, 0.56, 0.78])
            except Exception:
                pass
            ax.set_azimuth_ticks([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
            ax.set_azimuth_ticklabels(['0\u00b0', '30\u00b0', '60\u00b0', '90\u00b0',
                                        '120\u00b0', '150\u00b0', '180\u00b0', '210\u00b0',
                                        '240\u00b0', '270\u00b0', '300\u00b0', '330\u00b0'])
            # Match stereonet azimuth labels to the colourbar tick-label size
            # to keep spacing consistent around the net, especially when a
            # contour colour scale is displayed on the left.
            for _az_label in ax.get_azimuth_ticklabels():
                _az_label.set_fontsize(7)
            ax.grid(kind='equal_area_stereonet')

            pole_lines = None
            lin_lines = None
            artist_registry = defaultdict(list)
            stats_registry = {}
            label_registry = {}
            circle_registry = {}

            category_values = set()
            category_values.update(plane_categories)
            category_values.update(linear_categories)
            category_values.update(ref_categories)
            category_values.update([rec.get('category', 'All') for rec in kinematic_arrow_records])
            if not category_values:
                category_values.add('All')
            category_values = sorted(category_values, key=lambda value: str(value))
            category_styles = {
                category: self._default_category_style(i)
                for i, category in enumerate(category_values)
            }

            category_counts = {category: 0 for category in category_values}
            primary_categories = (linear_categories if show_linears and linear_categories
                                  else plane_categories if show_planes and plane_categories
                                  else ref_categories)
            for category in primary_categories:
                if category in category_counts:
                    category_counts[category] += 1

            # The classification legend and visibility controls are shown in
            # a right-hand Qt panel attached to the figure window. Keeping
            # these widgets outside the axes prevents overlap with the
            # stereonet and leaves the plot rendering unchanged.

            def _indices_for(categories, category):
                return [i for i, cat in enumerate(categories) if cat == category]

            for category in category_values:
                style = category_styles[category]

                if show_planes and has_planes:
                    idx = _indices_for(plane_categories, category)
                    if idx:
                        p_strikes = [plane_strikes[i] for i in idx]
                        p_dips = [plane_dips[i] for i in idx]
                        if stereoConfig.get('showGtCircles', False):
                            plane_artist = ax.plane(
                                p_strikes, p_dips, color=style.get('linecolor', style['color']),
                                linewidth=style['linewidth'], alpha=style['alpha'])
                            artist_registry[category].append({'artist': plane_artist, 'role': 'line'})
                        else:
                            pole_artist = ax.pole(
                                p_strikes, p_dips, linestyle='none',
                                marker=style['marker'], color=style['color'],
                                markersize=style['markersize'], alpha=style['alpha'])
                            artist_registry[category].append({'artist': pole_artist, 'role': 'marker'})

                        stats_by_mode, label_by_mode, circle_by_mode = _build_orientation_stat_artists(
                            ax, 'planes', p_strikes, p_dips, style)
                        if stats_by_mode is not None:
                            stats_registry[category] = stats_by_mode
                            label_registry[category] = label_by_mode
                            circle_registry[category] = circle_by_mode

                if show_linears and has_linears:
                    idx = _indices_for(linear_categories, category)
                    if idx:
                        l_plunges = [linear_plunges[i] for i in idx]
                        l_bearings = [linear_bearings[i] for i in idx]
                        line_artist = ax.line(
                            l_plunges, l_bearings, linestyle='none',
                            marker=style['marker'], color=style['color'],
                            markersize=style['markersize'], alpha=style['alpha'])
                        artist_registry[category].append({'artist': line_artist, 'role': 'marker'})

                        stats_by_mode, label_by_mode, circle_by_mode = _build_orientation_stat_artists(
                            ax, 'lines', l_plunges, l_bearings, style)
                        if stats_by_mode is not None:
                            stats_registry[category] = stats_by_mode
                            label_registry[category] = label_by_mode
                            circle_registry[category] = circle_by_mode

                if show_bearing_planes and effective_lin_planes:
                    idx = _indices_for(ref_categories, category)
                    if idx:
                        ref_artist = ax.plane(
                            [strikesref[i] for i in idx], [dipsref[i] for i in idx],
                            color=style.get('linecolor', style['color']), linewidth=style['linewidth'],
                            alpha=style['alpha'])
                        artist_registry[category].append({'artist': ref_artist, 'role': 'line'})
                    elif has_planes:
                        # Combined datasets such as Folds_PT may store the bearing
                        # plane in the regular Strike/Dip or DipDir/Dip fields.
                        idx = _indices_for(plane_categories, category)
                        if idx:
                            ref_artist = ax.plane(
                                [plane_strikes[i] for i in idx], [plane_dips[i] for i in idx],
                                color=style.get('linecolor', style['color']), linewidth=style['linewidth'],
                                alpha=style['alpha'])
                            artist_registry[category].append({'artist': ref_artist, 'role': 'line'})

                    if (plot_kinematics and kinematic_arrow_records and
                            effective_data_type == 'Lineations with Bearing Planes'):
                        arrow_records = [rec for rec in kinematic_arrow_records
                                         if rec.get('category', 'All') == category]
                        if arrow_records:
                            arrow_artists = self._plot_kinematic_arrows(
                                ax, arrow_records,
                                stereoConfig.get('kinematicsAnchor', 'Plane pole'),
                                color=style.get('arrowcolor', style['color']))
                            artist_registry[category].extend({'artist': arrow, 'role': 'arrow'} for arrow in (arrow_artists or []))

                if (plot_kinematics and kinematic_arrow_records and
                        not (show_bearing_planes and effective_lin_planes)):
                    arrow_records = [rec for rec in kinematic_arrow_records
                                     if rec.get('category', 'All') == category]
                    if arrow_records:
                        arrow_artists = self._plot_kinematic_arrows(
                            ax, arrow_records,
                            stereoConfig.get('kinematicsAnchor', 'Plane pole'),
                            color=style.get('arrowcolor', style['color']))
                        artist_registry[category].extend({'artist': arrow, 'role': 'arrow'} for arrow in (arrow_artists or []))

            contour_artists = []
            contour_colorbar = {'bar': None, 'cax': None, 'cid': None}

            def _position_contour_colorbar():
                """Keep the colour bar tied to the stereonet axes without resizing it."""
                cax = contour_colorbar.get('cax')
                if cax is None:
                    return
                pos = ax.get_position()
                width = 0.018
                # Keep the colourbar sufficiently far from stereonet azimuth
                # labels; this avoids crowding at ~240-300° during resizing.
                gap = 0.075
                height = pos.height * 0.72
                bottom = pos.y0 + (pos.height - height) / 2.0
                left = max(0.012, pos.x0 - gap - width)
                cax.set_position([left, bottom, width, height])

            def _remove_contour_colorbar():
                cid = contour_colorbar.get('cid')
                if cid is not None:
                    try:
                        ax.figure.canvas.mpl_disconnect(cid)
                    except Exception:
                        pass
                cbar = contour_colorbar.get('bar')
                if cbar is not None:
                    try:
                        cbar.remove()
                    except Exception:
                        pass
                cax = contour_colorbar.get('cax')
                if cax is not None:
                    try:
                        cax.remove()
                    except Exception:
                        pass
                contour_colorbar.update({'bar': None, 'cax': None, 'cid': None})

            def _add_contour_colorbar(contour_artist, measurement_label):
                """Add a continuous Kamb sigma-level colour scale.

                The colourbar is drawn in its own axes, so it never changes the
                stereonet axes size.  Its tick marks are the actual contour
                levels drawn by Matplotlib.
                """
                if contour_artist is None or contour_colorbar.get('bar') is not None:
                    return
                levels = np.asarray(getattr(contour_artist, 'levels', []), dtype=float)
                levels = levels[np.isfinite(levels)]
                if levels.size == 0:
                    return
                vmin, vmax = float(np.nanmin(levels)), float(np.nanmax(levels))
                if np.isclose(vmin, vmax):
                    vmax = vmin + 1.0
                cax = ax.figure.add_axes([0.02, 0.2, 0.018, 0.6])
                contour_colorbar['cax'] = cax
                _position_contour_colorbar()
                mappable = cm.ScalarMappable(
                    norm=Normalize(vmin=vmin, vmax=vmax),
                    cmap=getattr(contour_artist, 'cmap', cm.coolwarm))
                mappable.set_array([])
                cbar = ax.figure.colorbar(mappable, cax=cax)
                cbar.set_ticks(levels)
                cbar.set_ticklabels([f'{level:g}' for level in levels])

                # The colour bar is placed left of the stereonet.  Put both
                # tick labels and title on its exterior side so they do not
                # encroach on the stereonet circle.
                cbar.ax.yaxis.set_ticks_position('left')
                cbar.ax.yaxis.set_label_position('left')
                cbar.ax.tick_params(
                    labelsize=7, labelleft=True, labelright=False,
                    left=True, right=False, pad=3)
                cbar.set_label('Kamb contours (σ level)', fontsize=8, labelpad=8)
                contour_colorbar['bar'] = cbar
                contour_colorbar['cid'] = ax.figure.canvas.mpl_connect(
                    'resize_event', lambda event: _position_contour_colorbar())

            def _remove_contour_artist(contour_artist):
                if contour_artist is None:
                    return
                if hasattr(contour_artist, 'collections'):
                    for collection in list(contour_artist.collections):
                        try:
                            collection.remove()
                        except Exception:
                            collection.set_visible(False)
                elif isinstance(contour_artist, (list, tuple)):
                    for item in contour_artist:
                        _remove_contour_artist(item)
                elif hasattr(contour_artist, 'remove'):
                    try:
                        contour_artist.remove()
                    except Exception:
                        contour_artist.set_visible(False)

            def _update_visible_contours(visible_state=None):
                for contour_artist in list(contour_artists):
                    _remove_contour_artist(contour_artist)
                contour_artists[:] = []
                _remove_contour_colorbar()

                if not stereoConfig.get('showContours', True):
                    ax.figure.canvas.draw_idle()
                    return

                if visible_state is None:
                    visible_categories = set(category_values)
                else:
                    visible_categories = {category for category, visible in visible_state.items() if visible}

                if show_planes and has_planes:
                    idx = [i for i, category in enumerate(plane_categories)
                           if category in visible_categories]
                    if len(idx) >= 2:
                        contour_artist = ax.density_contour(
                            [plane_strikes[i] for i in idx],
                            [plane_dips[i] for i in idx],
                            measurement='poles', cmap=cm.coolwarm,
                            method='exponential_kamb', sigma=1.5,
                            linewidths=0.6, alpha=0.9, zorder=2)
                        contour_artists.append(contour_artist)
                        _add_contour_colorbar(contour_artist, 'Pole')

                if show_linears and has_linears:
                    idx = [i for i, category in enumerate(linear_categories)
                           if category in visible_categories]
                    if len(idx) >= 2:
                        contour_artist = ax.density_contour(
                            [linear_plunges[i] for i in idx],
                            [linear_bearings[i] for i in idx],
                            measurement='lines', cmap=cm.coolwarm,
                            method='exponential_kamb', sigma=1.5,
                            linewidths=0.6, alpha=0.9, zorder=2)
                        contour_artists.append(contour_artist)
                        _add_contour_colorbar(contour_artist, 'Lineation')

                ax.figure.canvas.draw_idle()

            girdle_artists = []

            def _update_visible_girdle(visible_state=None):
                for girdle_artist in list(girdle_artists):
                    _remove_contour_artist(girdle_artist)
                girdle_artists[:] = []

                if not stereoConfig.get('fitGirdle', False):
                    return
                if not (show_planes and has_planes):
                    return

                if visible_state is None:
                    visible_categories = set(category_values)
                else:
                    visible_categories = {category for category, visible in visible_state.items() if visible}

                idx = [i for i, category in enumerate(plane_categories)
                       if category in visible_categories]
                if len(idx) < 3:
                    return

                visible_strikes = [plane_strikes[i] for i in idx]
                visible_dips = [plane_dips[i] for i in idx]

                gs, gd = mplstereonet.fit_girdle(
                    visible_strikes, visible_dips, measurement='poles')
                girdle_artists.append(ax.plane(gs, gd, 'b-', linewidth=1.5))
                girdle_artists.append(ax.pole(gs, gd, 'ro', markersize=8))
                plunge, bearing = mplstereonet.pole2plunge_bearing(gs, gd)
                _, _, evals = mplstereonet.eigenvectors(
                    visible_strikes, visible_dips, measurement='poles')
                e1, e2, e3 = evals[0], evals[1], evals[2]
                pb = (f'Pole to best fit girdle\n'
                      f'  Plunge/Bearing: '
                      f'{int(round(plunge[0]))}/{int(round(bearing[0])):03d}')
                if e1 > e2 > e3 > 1e-10:
                    K = np.log(e1 / e2) / np.log(e2 / e3)
                    C = np.log(e1 / e3)
                    kshape = ('girdle' if K < 0.9 else
                              'cluster' if K > 1.1 else 'transitional')
                    info = (
                        f'{pb}\n'
                        f'  K = {K:.2f} ({kshape}),  C = {C:.2f}\n'
                        f'{"─" * 32}\n'
                        f'K (shape):\n'
                        f'  <1 = girdle  ≈1 = transitional  >1 = cluster\n'
                        f'C (strength):\n'
                        f'  0 = random → larger = stronger fabric'
                    )
                else:
                    info = pb
                girdle_artists.append(ax.text(
                    .9, -0.1, info, transform=ax.transAxes,
                    ha='left', va='bottom', fontsize=7, clip_on=False,
                    bbox=dict(boxstyle='round,pad=0.4', fc='lightblue',
                              ec='steelblue', alpha=0.85)))

            # Invisible reference artists used only by the existing interactive
            # selection logic. Visible category-specific artists are handled by
            # the category registry above.
            if show_linears and has_linears:
                lin_lines = ax.line(linear_plunges, linear_bearings,
                                    linestyle='none', marker='o', markersize=0,
                                    alpha=0)
            elif show_planes and has_planes and not stereoConfig.get('showGtCircles', False):
                pole_lines = ax.pole(plane_strikes, plane_dips,
                                     linestyle='none', marker='o', markersize=0,
                                     alpha=0)

            def _classification_legend_title():
                if not classification_field:
                    return 'Classification'
                for _layer in layers:
                    if _layer.type() != QgsMapLayer.VectorLayer:
                        continue
                    _idx = _layer.fields().lookupField(classification_field)
                    if _idx == -1:
                        continue
                    try:
                        _alias = _layer.attributeAlias(_idx) or ''
                    except Exception:
                        _alias = ''
                    return _alias if _alias else classification_field
                return classification_field or 'Classification'

            export_legend_artists = {}
            if classification_enabled and len(category_values) > 1:
                legend_handles = []
                legend_labels = []
                for category in category_values:
                    style = category_styles.get(category, self._default_category_style(0))
                    legend_handles.append(Line2D(
                        [0], [0], linestyle='none',
                        marker=style.get('marker', 'o'),
                        markerfacecolor=style.get('color', '#000000'),
                        markeredgecolor=style.get('color', '#000000'),
                        markersize=max(4.0, float(style.get('markersize', 5))),
                        alpha=float(style.get('alpha', 1.0))))
                    legend_labels.append(f'{category} (n={category_counts.get(category, 0)})')
                _ax_pos = ax.get_position()

                def _create_export_legend(ncol=1):
                    # Place the export legend to the right of the stereonet axes,
                    # not over the net.  Start with one column; only split into
                    # several columns if the rendered legend would extend below
                    # the bottom of the stereonet circle.
                    return fig.legend(
                        legend_handles, legend_labels,
                        title=_classification_legend_title(),
                        loc='upper left',
                        bbox_to_anchor=(min(0.985, _ax_pos.x1 + 0.025), _ax_pos.y1),
                        fontsize=7, title_fontsize=8, frameon=True,
                        borderpad=0.4, labelspacing=0.35, handletextpad=0.5,
                        columnspacing=0.8, ncol=max(1, int(ncol)))

                export_legend = _create_export_legend(1)

                try:
                    # Matplotlib can only know the true legend size after a draw.
                    # Use the rendered legend height in figure-fraction units and
                    # compare it with the stereonet axes height.  This keeps a
                    # single-column legend whenever it fits, and uses the minimum
                    # number of columns only when it would run below the net.
                    fig.canvas.draw()
                    renderer = fig.canvas.get_renderer()
                    legend_bbox = export_legend.get_window_extent(renderer=renderer)
                    fig_bbox = fig.bbox
                    legend_height_frac = legend_bbox.height / float(fig_bbox.height)
                    available_height_frac = max(0.05, _ax_pos.height)
                    if legend_height_frac > available_height_frac:
                        export_legend.remove()
                        legend_ncol = int(np.ceil(legend_height_frac / available_height_frac))
                        export_legend = _create_export_legend(legend_ncol)
                except Exception:
                    # Safe fallback for non-interactive backends.
                    pass
                export_legend.set_zorder(20)
                handles = getattr(export_legend, 'legend_handles', None)
                if handles is None:
                    handles = getattr(export_legend, 'legendHandles', [])
                texts = export_legend.get_texts()
                for category, handle, label_text in zip(category_values, handles, texts):
                    # Keep the legend label text black when styles are edited;
                    # only the marker handle should inherit the category colour.
                    label_text.set_color('black')
                    export_legend_artists[category] = [
                        {'artist': handle, 'role': 'marker'},
                        {'artist': label_text, 'role': 'legend_label'},
                    ]

            self._open_category_panel(
                fig, artist_registry,
                category_counts=category_counts,
                category_styles=category_styles,
                contour_update_callback=_update_visible_contours,
                girdle_update_callback=_update_visible_girdle,
                style_template_key=classification_field or 'All',
                export_legend_artists=export_legend_artists,
                show_visibility_controls=classification_enabled,
                stats_registry=stats_registry,
                label_registry=label_registry,
                circle_registry=circle_registry)

            if not classification_enabled:
                _visible_all = {category: True for category in category_values}
                _update_visible_contours(_visible_all)
                _update_visible_girdle(_visible_all)

            # Resolve which plotted points to use for interactive selection
            pts = None
            sel_fids = None
            sel_labels = None
            if lin_lines and linear_feature_ids:
                lon_data = np.array(lin_lines[0].get_xdata())
                lat_data = np.array(lin_lines[0].get_ydata())
                valid_pts = ~(np.isnan(lon_data) | np.isnan(lat_data))
                if valid_pts.any():
                    pts = np.column_stack([lon_data[valid_pts], lat_data[valid_pts]])
                    sel_fids   = [linear_feature_ids[i] for i, ok in enumerate(valid_pts) if ok]
                    sel_labels = [linear_labels[i]      for i, ok in enumerate(valid_pts) if ok]
            elif pole_lines and plane_feature_ids:
                lon_data = np.array(pole_lines[0].get_xdata())
                lat_data = np.array(pole_lines[0].get_ydata())
                valid_pts = ~(np.isnan(lon_data) | np.isnan(lat_data))
                if valid_pts.any():
                    pts = np.column_stack([lon_data[valid_pts], lat_data[valid_pts]])
                    sel_fids   = [plane_feature_ids[i] for i, ok in enumerate(valid_pts) if ok]
                    sel_labels = [plane_labels[i]      for i, ok in enumerate(valid_pts) if ok]

            if pts is not None and len(pts) > 0:
                sel_plot, = ax.plot([], [], 'ro', markersize=10, zorder=5,
                                    fillstyle='none', markeredgewidth=2)

                annot = ax.annotate(
                    '', xy=(0, 0), xycoords='data',
                    xytext=(8, 8), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow',
                              ec='gray', alpha=0.9),
                    fontsize=9, zorder=10)
                annot.set_visible(False)
                _hover_idx = [-1]

                def _on_hover(event):
                    if len(pts) == 0:
                        return
                    if event.inaxes != ax:
                        if _hover_idx[0] >= 0:
                            _hover_idx[0] = -1
                            annot.set_visible(False)
                            fig.canvas.draw_idle()
                        return
                    pts_axes = ax.transAxes.inverted().transform(
                        ax.transData.transform(pts))
                    cur_axes = ax.transAxes.inverted().transform(
                        [[event.x, event.y]])[0]
                    dists = np.hypot(pts_axes[:, 0] - cur_axes[0],
                                     pts_axes[:, 1] - cur_axes[1])
                    min_i = int(np.argmin(dists))
                    if dists[min_i] < 0.04:
                        if min_i != _hover_idx[0]:
                            _hover_idx[0] = min_i
                            annot.xy = (pts[min_i, 0], pts[min_i, 1])
                            annot.set_text(sel_labels[min_i])
                            annot.set_visible(True)
                            fig.canvas.draw_idle()
                    elif _hover_idx[0] >= 0:
                        _hover_idx[0] = -1
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

                _current_indices = [[]]
                _shift_held = [False]

                def _update_selection(indices):
                    _current_indices[0] = indices
                    if indices:
                        sel_plot.set_data(pts[indices, 0], pts[indices, 1])
                    else:
                        sel_plot.set_data([], [])
                    fig.canvas.draw_idle()
                    layer_sel = defaultdict(list)
                    for i in indices:
                        lyr, fid = sel_fids[i]
                        layer_sel[id(lyr)].append(fid)
                    all_layers = {id(lyr): lyr for lyr, _ in sel_fids}
                    for lid, lyr in all_layers.items():
                        lyr.selectByIds(layer_sel[lid]) if lid in layer_sel else lyr.removeSelection()

                def _on_lasso(verts):
                    pts_axes = ax.transAxes.inverted().transform(
                        ax.transData.transform(pts))
                    new_idx = np.where(Path(verts).contains_points(pts_axes))[0].tolist()
                    if _shift_held[0]:
                        combined = list(set(_current_indices[0]) | set(new_idx))
                    else:
                        combined = new_idx
                    _update_selection(combined)

                _press_xy = [None]

                def _on_press(event):
                    if event.inaxes == ax and event.button == 1:
                        _press_xy[0] = (event.x, event.y)

                def _on_release(event):
                    if event.button != 1 or _press_xy[0] is None or event.xdata is None:
                        _press_xy[0] = None
                        return
                    moved = (abs(event.x - _press_xy[0][0]) > 5 or
                             abs(event.y - _press_xy[0][1]) > 5)
                    _press_xy[0] = None
                    if moved or event.inaxes != ax:
                        return
                    dists = np.hypot(pts[:, 0] - event.xdata, pts[:, 1] - event.ydata)
                    min_i = int(np.argmin(dists))
                    if dists[min_i] < 0.1:
                        if event.key == 'shift':
                            combined = list(set(_current_indices[0]) | {min_i})
                        else:
                            combined = [min_i]
                        _update_selection(combined)

                def _on_key_press(event):
                    if event.key == 'escape':
                        _update_selection([])
                    elif event.key == 'shift':
                        _shift_held[0] = True

                def _on_key_release(event):
                    if event.key == 'shift':
                        _shift_held[0] = False

                self._stereonet_lasso = _BoundedLassoSelector(ax, _on_lasso)
                fig.canvas.mpl_connect('button_press_event', _on_press)
                fig.canvas.mpl_connect('button_release_event', _on_release)
                fig.canvas.mpl_connect('key_press_event', _on_key_press)
                fig.canvas.mpl_connect('key_release_event', _on_key_release)
                fig.canvas.mpl_connect('motion_notify_event', _on_hover)

                # Let the plot selection (lasso, click, or shift-click, all
                # wired above) double as a QGIS layer filter. Added directly
                # to the interactive category panel, which was already built
                # by _open_category_panel() above.
                #
                # This mirrors the "Filter Layer to Selected"/"Clear Filter"
                # routine from the Geochemistry Plotting Tools plugin: it
                # reads QGIS's own layer selection (already kept in sync by
                # _update_selection()'s selectByIds() calls above, so lasso/
                # click/shift-click selections apply here too) rather than
                # re-deriving it, and builds the subset string from the
                # primary key field's actual stored *value* per selected
                # feature - not the internal QGIS feature ID - which is what
                # PostGIS/Spatialite/GeoPackage layers require. Layers with no
                # declared primary key (e.g. Shapefile) fall back to OGR's
                # FID pseudo-column.
                controls = getattr(self, '_category_controls', None)
                panel = controls.get('panel') if controls else None
                if panel is not None and panel.layout() is not None:

                    def _build_fid_subset_string(lyr, feature_ids):
                        try:
                            pk_indexes = lyr.dataProvider().pkAttributeIndexes()
                        except Exception:
                            pk_indexes = []

                        if pk_indexes:
                            pk_field = lyr.fields()[pk_indexes[0]].name()
                            values = []
                            for fid in feature_ids:
                                feature = lyr.getFeature(fid)
                                if not feature.isValid():
                                    continue
                                val = feature[pk_field]
                                if isinstance(val, (int, float)):
                                    values.append(str(val))
                                else:
                                    values.append("'{}'".format(str(val).replace("'", "''")))
                            if values:
                                return '"{}" IN ({})'.format(pk_field, ','.join(values))

                        id_list = ','.join(str(fid) for fid in feature_ids)
                        return f"FID IN ({id_list})"

                    def _filter_to_selected():
                        all_layers = {id(lyr): lyr for lyr, _ in sel_fids}
                        filtered = []
                        failed = []
                        for lyr in all_layers.values():
                            selected_ids = lyr.selectedFeatureIds()
                            if not selected_ids:
                                continue
                            subset = _build_fid_subset_string(lyr, selected_ids)
                            if lyr.setSubsetString(subset):
                                # Re-select every (now sole) remaining feature.
                                # Some providers renumber feature IDs once a
                                # subset string is applied, which would
                                # otherwise silently drop or shift the
                                # selection - selectAll() sidesteps that, and
                                # guarantees the replot below picks up exactly
                                # the filtered features via selectedFeatures().
                                lyr.selectAll()
                                filtered.append(lyr.name())
                            else:
                                failed.append(lyr.name())
                        if not filtered and not failed:
                            QMessageBox.warning(
                                panel, 'Warning',
                                'No features are selected. Select points on the plot first '
                                '(lasso, click, or shift-click), then filter.')
                            return
                        if failed:
                            QMessageBox.warning(
                                panel, 'Filter failed',
                                "Could not filter this layer's data source to the selection "
                                f"for: {', '.join(failed)}.\nThis can happen with some data "
                                "providers. As an alternative, use QGIS's Export > Save Selected "
                                "Features As... to create a new layer from the selection.")
                        if filtered:
                            # Regenerate the plot from the now-filtered
                            # layer(s), so it shows only the selected/filtered
                            # entities instead of the stale, previously-plotted
                            # full set.
                            plt.close(fig)
                            self.contourPlot()

                    def _clear_filters():
                        all_layers = {id(lyr): lyr for lyr, _ in sel_fids}
                        cleared = [lyr.name() for lyr in all_layers.values()
                                   if lyr.subsetString() and lyr.setSubsetString('')]
                        if not cleared:
                            QMessageBox.information(
                                panel, 'No filter', "These layers aren't currently filtered.")

                    filter_group = QGroupBox('Layer Filter', panel)
                    filter_layout = QVBoxLayout(filter_group)
                    filter_caption = QLabel(
                        "Turn the current plot selection (lasso, click, or shift-click) "
                        "into a QGIS layer filter, restricting the layer(s) to just those "
                        "features.", filter_group)
                    filter_caption.setWordWrap(True)
                    filter_layout.addWidget(filter_caption)

                    filter_btn_row = QHBoxLayout()
                    filter_to_selected_btn = QPushButton('Filter Layer to Selected', filter_group)
                    filter_to_selected_btn.setToolTip(
                        "Set each affected layer's filter to only the currently selected features.")
                    clear_filter_btn = QPushButton('Clear Filter', filter_group)
                    clear_filter_btn.setToolTip('Remove the filter from every layer used in this plot.')
                    filter_btn_row.addWidget(filter_to_selected_btn)
                    filter_btn_row.addWidget(clear_filter_btn)
                    filter_layout.addLayout(filter_btn_row)

                    filter_to_selected_btn.clicked.connect(_filter_to_selected)
                    clear_filter_btn.clicked.connect(_clear_filters)

                    panel.layout().addWidget(filter_group)

            ax.set_title(layer.name() + " [# " + str(len(iter)) + "]", pad=24)
            plt.show()

        else:
            self.iface.messageBar().pushMessage("No data selected, or no structural data found: first select a layer with structural info, then select the points that you wish to plot", level=Qgis.Warning, duration=5)
        
    def fieldExists(self, layer, fieldnames):
        """Backward-compatible wrapper for older internal calls."""
        return self._field_exists(layer, fieldnames)
