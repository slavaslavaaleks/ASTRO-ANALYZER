# ASTRO-ANALYZER v0.1.0-beta.1
# Copyright (C) 2025 Алексеев Вячеслав Иванович
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import re
import time
import json
import math
import colorsys
import warnings
import traceback
import xml.etree.ElementTree as ET
from collections import deque

import queue
import threading
from PIL import Image, ImageTk
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import pandas as pd
import tifffile

import matplotlib.pyplot as plt

from scipy.ndimage import (
    gaussian_filter, binary_dilation, binary_opening,
    median_filter, convolve, uniform_filter, sobel,
    distance_transform_edt,
    label,
)

from skimage import measure
from skimage.transform import resize
from skimage.filters import threshold_otsu
from skimage.morphology import ball, remove_small_objects, binary_closing
from skimage.segmentation import watershed

warnings.filterwarnings('ignore')

TRANSLATIONS = {
    "ru": {
        "log.warn_prefix": "[WARN]",
        "log.loading": "    Загружаем: {name}",
        "log.2d_image": "    2D-изображение (Z=1) — будут 2D-параметры",
        "log.3d_stack": "    3D-стек: {n} срезов — будут 3D-параметры",
        "log.warn_only_one_channel": "[WARN] В файле только {n} канал(ов), DAPI задублирован",
        "log.warn_nd2_missing": "[WARN] nd2 не установлен, пробуем другие способы",
        "log.warn_nd2_failed": "[WARN] nd2 не смог прочитать: {err}",
        "log.warn_oiffile_missing": "[WARN] oiffile не установлен. Установите: pip install oiffile",
        "log.warn_oiffile_failed": "[WARN] oiffile не смог прочитать: {err}",
        "log.warn_oirfile_missing": "[WARN] oirfile не установлен. Установите: pip install oirfile",
        "log.warn_oirfile_failed": "[WARN] oirfile не смог прочитать: {err}",
        "log.warn_czifile_missing": "[WARN] czifile не установлен. Установите: pip install czifile",
        "log.warn_czifile_failed": "[WARN] czifile не смог прочитать: {err}",
        "log.warn_tifffile_failed": "[WARN] tifffile не смог прочитать: {err}",
        "log.warn_tifffile_imread_failed": "[WARN] tifffile.imread тоже не сработал: {err}",
        "log.warn_unknown_axes": "[WARN] Неизвестный порядок осей: {axes}",
        "log.warn_format_unsupported": "[WARN] Формат {ext} пока не поддерживается",
        "log.error": "Ошибка: {err}",
        "log.voxel_size_title": "    Размер вокселя",
        "log.file_label": "    Файл: {name}",
        "log.mode_manual": "    Режим: вручную",
        "log.scale_value": "    XY = {xy:.4f} µm/px, Z = {z:.4f} µm/срез",
        "log.no_voxel_metadata": "    Метаданные о размере пикселя отсутствуют",
        "log.use_panel_settings": "    Используются значения из настроек",
        "log.source_label": "    Источник: {src}",
        "log.xy_value": "    XY = {val:.4f} µm/px",
        "log.xy_not_defined": "    XY не определено",
        "log.z_value": "    Z = {val:.4f} µm/срез",
        "log.z_not_defined": "    Z не определено (возможно, 2D)",
        "log.channels_title": "\n    Каналы",
        "log.channel_line": "    канал {idx}, '{name}', эмиссия {em:.1f} нм, возбуждение {ex:.1f} нм",
        "log.channel_line_no_name": "    канал {idx}, эмиссия {em:.1f} нм, возбуждение {ex:.1f} нм",
        "log.file_read_failed": "    Не удалось прочитать файл",
        "log.channels_found": "    Найдено каналов: {n}",
        "log.single_channel": "    Один канал — используется как канал клеток, без DAPI",
        "log.asking_channel_choice": "    Запрашиваем выбор каналов у пользователя...",
        "log.channel_cancelled": "    Отмена — используются каналы по умолчанию",
        "log.cells_channel": "    Клетки: канал {idx}",
        "log.dapi_channel": "    DAPI: канал {idx}",
        "log.dapi_not_used": "    DAPI: не используется",
        "log.channels_in_analysis": "    Каналы в анализе: {ch}",
        "log.optics_title": "\n    Оптические параметры",
        "log.optics_manual": "    Режим: вручную",
        "log.optics_no_metadata": "    Оптических параметров в метаданных нет",
        "log.optics_cells": "    Клетки: NA = {na}, λ = {wl:.1f} нм",
        "log.optics_dapi": "    DAPI:   NA = {na}, λ = {wl:.1f} нм",
        "log.deconv_start": "\n    Выполняется 3D-деконволюция...",
        "log.deconv_dapi_start": "    Выполняется DAPI-деконволюция...",
        "log.deconv_slices_deleted": "    Удалено срезов по вашему выбору: {n} (Z: {z})",
        "log.deconv_all_slices_kept": "    Все срезы оставлены",
        "log.deconv_done": "    Деконволюция и ручной отбор завершены. Размер стэка: {shape}",
        "log.deconv_dapi_done": "    DAPI-деконволюция завершена. Размер стэка: {shape}",
        "log.soma_title": "\n    Разметка сом",
        "log.soma_no_points": "    Точки не добавлены. Сегментация отменена",
        "log.soma_points_set": "    Точек поставлено: {n}",
        "log.soma_cancelled": "    Расстановка отменена пользователем",
        "log.soma_starting": "    Сегментация сом...",
        "log.soma_done": "    Сегментация сом завершена",
        "log.nuclei_title": "    Сегментация ядер...",
        "log.nuclei_no_dapi": "    DAPI-канал не используется — ядра не сегментируются",
        "log.nuclei_no_soma": "    Сомы не найдены",
        "log.nuclei_found": "    Ядер найдено: {n}",
        "log.cells_title": "    Сегментация клеток...",
        "log.cells_no_soma": "    Сомы не найдены",
        "log.cells_no_signal": "    После вычитания фона сигнал отсутствует",
        "log.cells_no_growth": "    Областей для роста не найдено",
        "log.cells_found": "    Клеток найдено: {n}",
        "log.filter_title": "\n    Ручная фильтрация клеток",
        "log.filter_no_cells": "    Нет клеток для фильтрации",
        "log.filter_skipped": "    Фильтрация пропущена",
        "log.filter_no_selection": "    Проблемных клеток не выбрано",
        "log.filter_resegment_n": "    Пересегментация {n} клеток",
        "log.filter_review": "    Оценка результатов пересегментации",
        "log.filter_kept": "    Принято новых сегментов: {n}",
        "log.filter_reverted": "    Возвращено старых сегментов: {n}",
        "log.filter_deleted": "    Удалено клеток: {n}",
        "log.measure_2d_title": "\n    Измерение параметров (2D)",
        "log.measure_3d_title": "    Измерение параметров (3D)",
        "log.progress_done": "{prefix}Готово: {n}/{total} (100%) | всего {elapsed:.1f} с",
        "log.progress_update": "{prefix}Обработано {n}/{total} ({percent}%) | прошло {elapsed:.0f} с{eta}",
        "log.progress_eta": " | осталось {eta:.0f} с",
        "log.progress_closed_early": "{prefix}Завершено досрочно: {n}/{total} | всего {elapsed:.1f} с",
        "log.measurement_warn_shape": "    [WARN] {fn}: labels.shape={ls} != intensity.shape={is}. Возможно, downsample применён не везде.",
        "log.results_appended": "    Результаты добавлены в {path} (всего строк: {n})",
        "log.results_read_failed": "    Не удалось прочитать существующий {path}: {err}",
        "log.viz_overlay_saved": "    Overlay-стек сохранён: {name}",
        "log.viz_saved_n_cells": "    Сохранено {n} клеток",
        "log.viz_final_start": "\n    Генерация финальных визуализаций...",
        "log.viz_saved_to": "    Визуализации сохранены в {path}",
        "log.viz_deconv_cells_saved": "    Сохранено деконволюированное изображение КЛЕТОК: {shape}",
        "log.viz_deconv_dapi_saved": "    Сохранено деконволюированное DAPI: {shape}",
        "log.viz_mip_cells_failed": "    Не удалось сохранить MIP КЛЕТОК: {err}",
        "log.viz_mip_dapi_failed": "    Не удалось сохранить MIP DAPI: {err}",
        "log.slice_title": "Z-срез: {z}",
        "log.deconv_cells_title": "{name} — КЛЕТКИ (деконволюция, MIP)",
        "log.deconv_dapi_title": "{name} — DAPI (деконволюция, MIP)",
        "log.title_soma_slices": "Сегментация сом",
        "log.title_nuclei_slices": "Сегментация ядер",
        "log.title_cells_slices": "Визуализация целой клетки с отростками",
        "log.batch_type_cells": "    Тип клеток: {t}",
        "log.batch_image_mode": "    Тип изображения: {m}",
        "log.batch_deconv": "    Деконволюция: {s}",
        "log.batch_filter": "    Ручная фильтрация: {s}",
        "log.batch_files_count": "    Файлов к обработке: {n}",
        "log.batch_file_n": "\n    Обработка файла {i} из {total}: {name}",
        "log.batch_files_done": "    Файлов обработано: {n}",
        "log.batch_cells_found": "    Клеток найдено: {n}",
        "log.batch_total_time": "    Общее время: {t:.1f} с",
        "log.batch_averages": "\n    Средние значения по сессии",
        "log.batch_avg_area": "    Площадь клетки: {v:.0f} мкм²",
        "log.batch_avg_nuc_area": "    Площадь ядра: {v:.0f} мкм²",
        "log.batch_avg_volume": "    Объём клетки: {v:.0f} мкм³",
        "log.batch_avg_nuc_volume": "    Объём ядра: {v:.0f} мкм³",
        "log.batch_avg_nc_ratio": "    Ядро/цитоплазма: {v:.2f}",
        "log.batch_avg_elong": "    Элонгация: {v:.2f}",
        "log.batch_no_cells": "\n    Клетки не найдены ни в одном файле",
        "log.file_summary": "\n    Итог по файлу",
        "log.file_summary_cells": "    Клеток найдено: {n}",
        "log.file_summary_time": "    Общее время: {t:.1f} с",
        "log.dapi_not_used_msg": "    DAPI не используется — ядра сегментироваться не будут",
        "log.soma_not_marked": "    Сомы не размечены — файл пропущен",
        "log.filter_skip_no_norm": "    Ручная фильтрация пропущена (нет astro_norm)",
        "log.filter_skip": "    Ручная фильтрация пропущена",
        "log.segmentation_time": "    Время сегментации: {t:.1f} с",
        "log.measurement_time": "    Время измерения: {t:.1f} с",
        "log.error_prefix": "Ошибка: {err}",
        "log.image_size": "    Размер изображения: {w}×{h} → даунсемплинг в {n} раза",
        "log.image_size_failed": "    Не удалось определить размер изображения: {err}",
        "log.scale_from_panel": "    Масштаб из настроек: XY = {xy:.4f} µm, Z = {z:.4f} µm",
        "log.analysis_done": "    Анализ завершён",
        "log.analysis_done_title": "Анализ завершён",
        "log.metadata_reading": "    Чтение метаданных: {name}",
        "log.metadata_xy_ok": "    XY = {v:.4f} µm/px",
        "log.metadata_xy_none": "    XY не найдено — оставлено значение из настроек",
        "log.metadata_z_ok": "    Z = {v:.4f} µm/срез",
        "log.metadata_z_none": "    Z не найдено — оставлено значение из настроек",
        "log.metadata_fields_updated": "    Поля настроек обновлены",
        "log.metadata_optics_read_error": "Ошибка чтения оптики: {err}",
        "log.metadata_scale_read_error": "Ошибка чтения масштаба: {err}",

        "cell.astrocytes": "астроцитов",
        "cell.microglia": "микроглии",
        "cell.astrocytes_label": "Астроциты",
        "cell.microglia_label": "Микроглия",
        "cell.neurons": "нейронов",
        "cell.neurons_label": "Нейроны",
        "mode.confocal": "Конфокальный",
        "mode.epifluorescence": "Эпифлуоресцентный",
        "mode.confocal_full": "конфокальный",
        "mode.epifluorescence_full": "эпифлуоресцентный",
        "state.on": "включена",
        "state.off": "выключена",

        "ui.settings.title_main": "Укажите файл и параметры — затем «Запустить анализ».",
        "ui.settings.subtitle": "После этого будет проведена сегментация и морфометрия клеток на фотографии.",
        "ui.settings.file": "Файл/папка:",
        "ui.settings.browse": "Обзор…",
        "ui.settings.load_mode": "Режим загрузки:",
        "ui.settings.load_one": "Один файл",
        "ui.settings.load_dir": "Папка",
        "ui.settings.metadata": "Метаданные:",
        "ui.settings.read_metadata": "Прочитать из файла",
        "ui.settings.cell_type": "Тип клеток:",
        "ui.settings.image_mode": "Тип изображения:",
        "ui.settings.deconv": "Деконволюция:",
        "ui.settings.use_deconv": "Применять 3D-деконволюцию перед анализом",
        "ui.settings.xy": "XY (мкм/пиксель):",
        "ui.settings.z": "Z (мкм/срез):",
        "ui.settings.scale_hint": "Введите масштаб сюда, либо нажмите «Прочитать из файла» — значения подставятся сами.",
        "ui.settings.cells": "Клетки:",
        "ui.settings.dapi": "DAPI:",
        "ui.settings.default": "[по умолчанию]",
        "ui.settings.from_metadata": "[из метаданных]",
        "ui.settings.filter": "Ручная фильтрация клеток:",
        "ui.settings.run_filter": "Запускать после сегментации",
        "ui.settings.run": "▶   Запустить анализ",
        "ui.settings.exit": "Выход",
        "ui.settings.back_to_log": "→ Вернуться к логу",
        "ui.settings.back_to_settings": "← Настройки",
        "ui.settings.log_title": "Протокол работы",

        "ui.header.lang": "Язык:",

        "ui.settings.na": "NA:",
        "ui.settings.lambda_nm": "λ(нм):",
        "ui.settings.epsilon": "ε:",

        "btn.ok": "OK",
        "btn.cancel": "Отмена",
        "btn.apply": "Применить",
        "btn.done": "Готово",
        "btn.continue": "Продолжить",
        "btn.skip": "Пропустить",
        "btn.start_analysis": "Запустить сегментацию",
        "btn.choose": "Выбрать…",
        "btn.delete_return_slice": "Удалить / Вернуть срез",

        "dlg.files.title": "Выбор файлов",
        "dlg.files.one_or_dir": "Один файл или папка?",
        "dlg.files.one": "Один файл",
        "dlg.files.dir": "Папка с файлами",
        "dlg.files.selected_n": "Выбрано: {n} файл(ов)",
        "dlg.files.nothing": "Ничего не выбрано.",
        "dlg.files.images": "Изображения",
        "dlg.files.all": "Все файлы",

        "dlg.channels.title": "Выбор каналов",
        "dlg.channels.instructions": "Отметьте каналы, которые войдут в анализ",
        "dlg.channels.no_dapi_checkbox": "Нет DAPI-канала (ядра не сегментировать)",
        "dlg.channels.channel_n": "КАНАЛ {n}",
        "dlg.channels.role_cells": "КЛЕТКИ",
        "dlg.channels.role_dapi": "DAPI",
        "dlg.channels.role_none": "—",
        "dlg.channels.take": "✓ брать в анализ",
        "dlg.channels.skip": "не брать",
        "dlg.channels.choose_one": "Выберите хотя бы один канал.",
        "dlg.channels.help": "ЛКМ по картинке — назначить канал КЛЕТОК.\nПКМ по картинке — назначить канал DAPI.\nЛКМ по нижней кнопке — вкл/выкл канал.",

        "dlg.dapi.title": "Определение DAPI-канала",
        "dlg.dapi.text": "Автоматически определить не удалось.\nКликните по каналу, который является DAPI:",
        "dlg.dapi.help": "Кликните по каналу с ядрами.",

        "dlg.soma.title": "Расстановка сом",
        "dlg.soma.counter": "Точек: {n}",
        "dlg.soma.coords": "XY = ({y}, {x})",
        "dlg.soma.help": "ЛКМ — добавить сому,\nCtrl+ЛКМ — удалить ближайшую,\nПКМ — открыть лупу (оверлей).\nВ оверлее: слайдер Z, ЛКМ — поставить сому.",

        "dlg.soma.overlay_z": "Z = {z}",

        "dlg.soma_confirm.title": "Подтверждение расстановки сом",
        "dlg.soma_confirm.z_slice": "Z-срез",
        "dlg.soma_confirm.z_badge": "Z = {z}   ·   сом здесь: {n}",
        "dlg.soma_confirm.info": "Всего сом: {total}   ·   на этом срезе: {here}",
        "dlg.soma_confirm.help": "Листайте по всем Z.\nЗелёные точки — центры сом.",

        "dlg.slices.title": "Отбор Z-срезов",
        "dlg.slices.z_slice": "Z-срез",
        "dlg.slices.kept": "ОСТАВЛЕН",
        "dlg.slices.deleted": "УДАЛЁН",
        "dlg.slices.badge": "Z = {z}   [{label}]",
        "dlg.slices.toggle": "Удалить / Вернуть срез",
        "dlg.slices.help": "Листайте слайдером.\nКнопка переключает текущий срез.",
        "dlg.slices.title_selector": "Кликните на срезы для удаления (X). Enter — применить.",

        "dlg.filter.title": "Ручная фильтрация клеток",
        "dlg.filter.help": "ЛКМ — НЕДОРОСТ (синий),\nПКМ — ПЕРЕРОСТ (красный).",

        "dlg.review.title": "Ревью пересегментации",
        "dlg.review.help": "ЛКМ — ПРИНЯТЬ (зелёный),\nПКМ — УДАЛИТЬ (красный).",

        "dlg.repick.title": "Повторная разметка сомы",
        "dlg.repick.point_warning": "Поставьте точку ЛКМ.",
        "dlg.repick.help": "Красный контур — старый сегмент.\nЛКМ — поставить точку.",

        "dlg.unknown": "Неизвестный диалог: {kind}",
        "log.dialog_error": "ОШИБКА ДИАЛОГА: {err}",

        "dlg.error.title": "Ошибка",
        "dlg.files.title_short": "Файлы",
        "dlg.files.no_supported": "Поддерживаемые файлы не найдены.",

        "app.title": "ASTRO-ANALYZER",
        "app.header_title": "3D анализ клеток нервной ткани",
    },

    "en": {
        "log.warn_prefix": "[WARN]",
        "log.loading": "    Loading: {name}",
        "log.2d_image": "    2D image (Z=1) — measuring 2D parameters",
        "log.3d_stack": "    3D stack: {n} slices — measuring 3D parameters",
        "log.warn_only_one_channel": "[WARN] File contains only {n} channel(s), DAPI duplicated",
        "log.warn_nd2_missing": "[WARN] nd2 not installed, trying other readers",
        "log.warn_nd2_failed": "[WARN] nd2 failed: {err}",
        "log.warn_czifile_missing": "[WARN] czifile not installed. Install: pip install czifile",
        "log.warn_czifile_failed": "[WARN] czifile failed: {err}",
        "log.warn_oiffile_missing": "[WARN] oiffile not installed. Install: pip install oiffile",
        "log.warn_oiffile_failed": "[WARN] oiffile failed: {err}",
        "log.warn_oirfile_missing": "[WARN] oirfile not installed. Install: pip install oirfile",
        "log.warn_oirfile_failed": "[WARN] oirfile failed: {err}",
        "log.warn_tifffile_failed": "[WARN] tifffile failed: {err}",
        "log.warn_tifffile_imread_failed": "[WARN] tifffile.imread also failed: {err}",
        "log.warn_unknown_axes": "[WARN] Unknown axis order: {axes}",
        "log.warn_format_unsupported": "[WARN] Format {ext} is not supported yet",
        "log.error": "Error: {err}",
        "log.voxel_size_title": "    Voxel size",
        "log.file_label": "    File: {name}",
        "log.mode_manual": "    Mode: manual",
        "log.scale_value": "    XY = {xy:.4f} µm/px, Z = {z:.4f} µm/slice",
        "log.no_voxel_metadata": "    No pixel-size metadata found",
        "log.use_panel_settings": "    Using values from the settings panel",
        "log.source_label": "    Source: {src}",
        "log.xy_value": "    XY = {val:.4f} µm/px",
        "log.xy_not_defined": "    XY not defined",
        "log.z_value": "    Z = {val:.4f} µm/slice",
        "log.z_not_defined": "    Z not defined (possibly 2D)",
        "log.channels_title": "\n    Channels",
        "log.channel_line": "    channel {idx}, '{name}', emission {em:.1f} nm, excitation {ex:.1f} nm",
        "log.channel_line_no_name": "    channel {idx}, emission {em:.1f} nm, excitation {ex:.1f} nm",
        "log.file_read_failed": "    Failed to read the file",
        "log.channels_found": "    Channels found: {n}",
        "log.single_channel": "    Single channel — used as cell channel, no DAPI",
        "log.asking_channel_choice": "    Asking user to choose channels...",
        "log.channel_cancelled": "    Cancelled — using default channels",
        "log.cells_channel": "    Cells: channel {idx}",
        "log.dapi_channel": "    DAPI: channel {idx}",
        "log.dapi_not_used": "    DAPI: not used",
        "log.channels_in_analysis": "    Channels in analysis: {ch}",
        "log.optics_title": "\n    Optical parameters",
        "log.optics_manual": "    Mode: manual",
        "log.optics_no_metadata": "    No optical parameters in metadata",
        "log.optics_cells": "    Cells: NA = {na}, λ = {wl:.1f} nm",
        "log.optics_dapi": "    DAPI:  NA = {na}, λ = {wl:.1f} nm",
        "log.deconv_start": "\n    Running 3D deconvolution...",
        "log.deconv_dapi_start": "    Running DAPI deconvolution...",
        "log.deconv_slices_deleted": "    Slices removed by your choice: {n} (Z: {z})",
        "log.deconv_all_slices_kept": "    All slices kept",
        "log.deconv_done": "    Deconvolution and manual selection complete. Stack shape: {shape}",
        "log.deconv_dapi_done": "    DAPI deconvolution complete. Stack shape: {shape}",
        "log.soma_title": "\n    Soma marking",
        "log.soma_no_points": "    No points added. Segmentation cancelled",
        "log.soma_points_set": "    Points set: {n}",
        "log.soma_cancelled": "    Marking cancelled by user",
        "log.soma_starting": "    Segmenting somas...",
        "log.soma_done": "    Soma segmentation complete",
        "log.nuclei_title": "    Segmenting nuclei...",
        "log.nuclei_no_dapi": "    DAPI channel not used — nuclei are not segmented",
        "log.nuclei_no_soma": "    No somas found",
        "log.nuclei_found": "    Nuclei found: {n}",
        "log.cells_title": "    Segmenting cells...",
        "log.cells_no_soma": "    No somas found",
        "log.cells_no_signal": "    No signal left after background subtraction",
        "log.cells_no_growth": "    No regions available for growth",
        "log.cells_found": "    Cells found: {n}",
        "log.filter_title": "\n    Manual cell filtering",
        "log.filter_no_cells": "    No cells to filter",
        "log.filter_skipped": "    Filtering skipped",
        "log.filter_no_selection": "    No problematic cells selected",
        "log.filter_resegment_n": "    Re-segmenting {n} cells",
        "log.filter_review": "    Reviewing re-segmentation results",
        "log.filter_kept": "    New segments kept: {n}",
        "log.filter_reverted": "    Old segments reverted: {n}",
        "log.filter_deleted": "    Cells deleted: {n}",
        "log.measure_2d_title": "\n    Measuring parameters (2D)",
        "log.measure_3d_title": "    Measuring parameters (3D)",
        "log.progress_done": "{prefix}Done: {n}/{total} (100%) | total {elapsed:.1f} s",
        "log.progress_update": "{prefix}Processed {n}/{total} ({percent}%) | elapsed {elapsed:.0f} s{eta}",
        "log.progress_eta": " | ETA {eta:.0f} s",
        "log.progress_closed_early": "{prefix}Closed early: {n}/{total} | total {elapsed:.1f} s",
        "log.measurement_warn_shape": "    [WARN] {fn}: labels.shape={ls} != intensity.shape={is}. Downsample may not be applied everywhere.",
        "log.results_appended": "    Results appended to {path} (total rows: {n})",
        "log.results_read_failed": "    Failed to read existing {path}: {err}",
        "log.viz_overlay_saved": "    Overlay stack saved: {name}",
        "log.viz_saved_n_cells": "    Saved {n} cells",
        "log.viz_final_start": "\n    Generating final visualizations...",
        "log.viz_saved_to": "    Visualizations saved to {path}",
        "log.viz_deconv_cells_saved": "    Deconvolved CELLS image saved: {shape}",
        "log.viz_deconv_dapi_saved": "    Deconvolved DAPI saved: {shape}",
        "log.viz_mip_cells_failed": "    Failed to save CELLS MIP: {err}",
        "log.viz_mip_dapi_failed": "    Failed to save DAPI MIP: {err}",
        "log.slice_title": "Z-slice: {z}",
        "log.deconv_cells_title": "{name} — CELLS (deconvolution, MIP)",
        "log.deconv_dapi_title": "{name} — DAPI (deconvolution, MIP)",
        "log.title_soma_slices": "Soma segmentation",
        "log.title_nuclei_slices": "Nuclei segmentation",
        "log.title_cells_slices": "Whole cell visualization with processes",
        "log.batch_type_cells": "    Cell type: {t}",
        "log.batch_image_mode": "    Image type: {m}",
        "log.batch_deconv": "    Deconvolution: {s}",
        "log.batch_filter": "    Manual filtering: {s}",
        "log.batch_files_count": "    Files to process: {n}",
        "log.batch_file_n": "\n    Processing file {i} of {total}: {name}",
        "log.batch_files_done": "    Files processed: {n}",
        "log.batch_cells_found": "    Cells found: {n}",
        "log.batch_total_time": "    Total time: {t:.1f} s",
        "log.batch_averages": "\n    Session averages",
        "log.batch_avg_area": "    Cell area: {v:.0f} µm²",
        "log.batch_avg_nuc_area": "    Nucleus area: {v:.0f} µm²",
        "log.batch_avg_volume": "    Cell volume: {v:.0f} µm³",
        "log.batch_avg_nuc_volume": "    Nucleus volume: {v:.0f} µm³",
        "log.batch_avg_nc_ratio": "    Nucleus/cytoplasm: {v:.2f}",
        "log.batch_avg_elong": "    Elongation: {v:.2f}",
        "log.batch_no_cells": "\n    No cells found in any file",
        "log.file_summary": "\n    File summary",
        "log.file_summary_cells": "    Cells found: {n}",
        "log.file_summary_time": "    Total time: {t:.1f} s",
        "log.dapi_not_used_msg": "    DAPI not used — nuclei will not be segmented",
        "log.soma_not_marked": "    Somas not marked — file skipped",
        "log.filter_skip_no_norm": "    Manual filtering skipped (no astro_norm)",
        "log.filter_skip": "    Manual filtering skipped",
        "log.segmentation_time": "    Segmentation time: {t:.1f} s",
        "log.measurement_time": "    Measurement time: {t:.1f} s",
        "log.error_prefix": "Error: {err}",
        "log.image_size": "    Image size: {w}×{h} → downsampling by {n}",
        "log.image_size_failed": "    Failed to determine image size: {err}",
        "log.scale_from_panel": "    Scale from settings: XY = {xy:.4f} µm, Z = {z:.4f} µm",
        "log.analysis_done": "    Analysis complete",
        "log.analysis_done_title": "Analysis complete",
        "log.metadata_reading": "    Reading metadata: {name}",
        "log.metadata_xy_ok": "    XY = {v:.4f} µm/px",
        "log.metadata_xy_none": "    XY not found — keeping value from settings",
        "log.metadata_z_ok": "    Z = {v:.4f} µm/slice",
        "log.metadata_z_none": "    Z not found — keeping value from settings",
        "log.metadata_fields_updated": "    Settings fields updated",
        "log.metadata_optics_read_error": "Error reading optics: {err}",
        "log.metadata_scale_read_error": "Error reading scale: {err}",

        "cell.astrocytes": "astrocytes",
        "cell.microglia": "microglia",
        "cell.astrocytes_label": "Astrocytes",
        "cell.microglia_label": "Microglia",
        "cell.neurons": "neurons",
        "cell.neurons_label": "Neurons",
        "mode.confocal": "Confocal",
        "mode.epifluorescence": "Epifluorescence",
        "mode.confocal_full": "confocal",
        "mode.epifluorescence_full": "epifluorescence",
        "state.on": "enabled",
        "state.off": "disabled",

        "ui.settings.title_main": "Choose a file and parameters, then click \"Run analysis\".",
        "ui.settings.subtitle": "After that, cell segmentation and morphometry will be performed.",
        "ui.settings.file": "File/folder:",
        "ui.settings.browse": "Browse…",
        "ui.settings.load_mode": "Load mode:",
        "ui.settings.load_one": "Single file",
        "ui.settings.load_dir": "Folder",
        "ui.settings.metadata": "Metadata:",
        "ui.settings.read_metadata": "Read from file",
        "ui.settings.cell_type": "Cell type:",
        "ui.settings.image_mode": "Image type:",
        "ui.settings.deconv": "Deconvolution:",
        "ui.settings.use_deconv": "Apply 3D deconvolution before analysis",
        "ui.settings.xy": "XY (µm/pixel):",
        "ui.settings.z": "Z (µm/slice):",
        "ui.settings.scale_hint": "Enter the scale here, or click \"Read from file\" — values will be filled automatically.",
        "ui.settings.cells": "Cells:",
        "ui.settings.dapi": "DAPI:",
        "ui.settings.default": "[default]",
        "ui.settings.from_metadata": "[from metadata]",
        "ui.settings.filter": "Manual cell filtering:",
        "ui.settings.run_filter": "Run after segmentation",
        "ui.settings.run": "▶   Run analysis",
        "ui.settings.exit": "Exit",
        "ui.settings.back_to_log": "→ Back to log",
        "ui.settings.back_to_settings": "← Settings",
        "ui.settings.log_title": "Log",

        "ui.header.lang": "Lang:",

        "ui.settings.na": "NA:",
        "ui.settings.lambda_nm": "λ(nm):",
        "ui.settings.epsilon": "ε:",

        "btn.ok": "OK",
        "btn.cancel": "Cancel",
        "btn.apply": "Apply",
        "btn.done": "Done",
        "btn.continue": "Continue",
        "btn.skip": "Skip",
        "btn.start_analysis": "Start segmentation",
        "btn.choose": "Choose…",
        "btn.delete_return_slice": "Delete / Restore slice",

        "dlg.files.title": "File selection",
        "dlg.files.one_or_dir": "Single file or folder?",
        "dlg.files.one": "Single file",
        "dlg.files.dir": "Folder with files",
        "dlg.files.selected_n": "Selected: {n} file(s)",
        "dlg.files.nothing": "Nothing selected.",
        "dlg.files.images": "Images",
        "dlg.files.all": "All files",

        "dlg.channels.title": "Channel selection",
        "dlg.channels.instructions": "Mark the channels to include in the analysis",
        "dlg.channels.no_dapi_checkbox": "No DAPI channel (do not segment nuclei)",
        "dlg.channels.channel_n": "CHANNEL {n}",
        "dlg.channels.role_cells": "CELLS",
        "dlg.channels.role_dapi": "DAPI",
        "dlg.channels.role_none": "—",
        "dlg.channels.take": "✓ include",
        "dlg.channels.skip": "skip",
        "dlg.channels.choose_one": "Select at least one channel.",
        "dlg.channels.help": "Left click image — assign CELLS channel.\nRight click image — assign DAPI channel.\nLeft click bottom button — toggle channel.",

        "dlg.dapi.title": "DAPI channel detection",
        "dlg.dapi.text": "Automatic detection failed.\nClick the channel that contains DAPI:",
        "dlg.dapi.help": "Click the channel with nuclei.",

        "dlg.soma.title": "Soma marking",
        "dlg.soma.counter": "Points: {n}",
        "dlg.soma.coords": "XY = ({y}, {x})",
        "dlg.soma.help": "LMB — add soma,\nCtrl+LMB — delete nearest,\nRMB — open magnifier (overlay).\nIn overlay: Z slider, LMB — place soma.",

        "dlg.soma.overlay_z": "Z = {z}",

        "dlg.soma_confirm.title": "Confirm soma marking",
        "dlg.soma_confirm.z_slice": "Z-slice",
        "dlg.soma_confirm.z_badge": "Z = {z}   ·   somas here: {n}",
        "dlg.soma_confirm.info": "Total somas: {total}   ·   on this slice: {here}",
        "dlg.soma_confirm.help": "Browse all Z.\nGreen dots — soma centers.",

        "dlg.slices.title": "Z-slice selection",
        "dlg.slices.z_slice": "Z-slice",
        "dlg.slices.kept": "KEPT",
        "dlg.slices.deleted": "DELETED",
        "dlg.slices.badge": "Z = {z}   [{label}]",
        "dlg.slices.toggle": "Delete / Restore slice",
        "dlg.slices.help": "Use the slider.\nThe button toggles the current slice.",
        "dlg.slices.title_selector": "Click on slices to delete (X). Enter — apply.",

        "dlg.filter.title": "Manual cell filtering",
        "dlg.filter.help": "LMB — UNDERSIZED (blue),\nRMB — OVERSIZED (red).",

        "dlg.review.title": "Re-segmentation review",
        "dlg.review.help": "LMB — ACCEPT (green),\nRMB — DELETE (red).",

        "dlg.repick.title": "Re-mark soma",
        "dlg.repick.point_warning": "Place a point with LMB.",
        "dlg.repick.help": "Red outline — old segment.\nLMB — place a point.",

        "dlg.unknown": "Unknown dialog: {kind}",
        "log.dialog_error": "DIALOG ERROR: {err}",

        "dlg.error.title": "Error",
        "dlg.files.title_short": "Files",
        "dlg.files.no_supported": "No supported files found.",

        "app.title": "ASTRO-ANALYZER",
        "app.header_title": "3D analysis of nerve tissue cells",
    }
}

class AnalysisContext:
    def __init__(self):
        self.astrocyte_channel = 0
        self.dapi_channel = 1
        self.cell_type = "astrocytes"
        self.cell_type_ru = None
        self.is_2d = False
        self.bridge = None

        self.lang = "ru"

        self.image_mode = "confocal"

        self._last_astro_norm = None

        self.soma_segmented = None
        self._soma_labels = None
        self.nuclei_segmented = None
        self.cells_segmented = None

        self._channels_initialized = False
        self._optics_initialized = False

        self.downsample = 1
        self.growth_threshold = 1.3
        self.background_percentile = 15
        self.min_nucleus_size = 10
        self.max_growth_z = 4
        self.dilation_iterations = 5
        self.min_cell_volume = 100
        self.surface_ratio_threshold = 5.0
        self.elongation_threshold = 2.0
        self.solidity_threshold = 0.9
        self.intensity_ratio = 0.8
        self.min_periphery_size = 50
        self.nucleus_smoothing = 2.0
        self.nucleus_min_distance = 8
        self.nucleus_threshold_factor = 0.4
        self.nucleus_min_size = 50
        self.detail_min_process_length = 5
        self.detail_min_branch_length = 3
        self.detail_intensity_percentiles = [10, 25, 50, 75, 90]
        self.detail_distance_bins = 10
        self.auto_mode = True
        self.adaptive_params = {
            'roughness_threshold': 1.15,
            'circularity_threshold': 1.5,
            'solidity_threshold': 0.92,
            'sharp_ratio_threshold': 0.25,
            'filter_strictness': 1.0
        }
        self.xy_scale_um = None
        self.z_scale_um = None
        self.scale_initialized = False
        self.nucleus_seed_iterations = 3
        self.min_seed_mass_ratio = 1.8
        self.min_seed_cell_volume_ratio = 4.0
        self.seed_gfap_percentile = 60
        self.seed_gfap_min_fraction = 0.85
        self.seed_local_shell_dilation = 3
        self.deconvolve_images = True
        self._deconvolution_applied = False
        self._deconv_params = {
            'na': 0.95,
            'emission_wavelength_um': 0.6995,
            'wiener_epsilon': 0.005
        }
        self._dapi_deconv_params = {
            'na': 0.95,
            'emission_wavelength_um': 0.460,
            'wiener_epsilon': 0.005
        }
        self._forced_dapi_idx = None
        self._forced_astro_idx = None

        self._current_base_name = None
        self._current_out_dir = None

        self._gui_channel_mode = "auto"
        self._gui_scale_mode = "auto"
        self._gui_deconv_mode = "auto"
        self._gui_do_filter = True
        self._gui_settings_locked = False
        self._gui_metadata_cache = {}
        self._gui_scale_confirmed = False
        self._gui_optics_confirmed = False

        self.has_dapi = True
        self.n_channels_original = 1
        self.cell_channel_idx = None
        self.dapi_channel_idx = None
        self.selected_channels = None

        self._input_xy_um = None
        self._input_z_um = None

    def t(self, key, **kwargs):
        lang = getattr(self, "lang", "ru")
        table = TRANSLATIONS.get(lang, TRANSLATIONS["ru"])
        template = table.get(key)
        if template is None:
            template = TRANSLATIONS["ru"].get(key, key)
        try:
            return template.format(**kwargs)
        except Exception:
            return template

    def log(self, *args, **kwargs):
        msg = " ".join(str(a) for a in args)
        if self.bridge is not None:
            self.bridge.log(msg)
        else:
            print(msg, **kwargs)

def compute_downsample(y_dim, x_dim, target=512):
    bigger = max(int(y_dim), int(x_dim))
    if bigger <= target:
        return 1
    return int(math.ceil(bigger / target))

class FileIO:
    def __init__(self, ctx):
        self.ctx = ctx

    def _read_image_universal(self, filepath):
        ext = os.path.splitext(filepath)[1].lower()

        if ext in ['.oib', '.oif']:
            try:
                import oiffile
                with oiffile.OifFile(filepath) as oib:
                    arr = oib.asarray()
                    axes = getattr(oib, 'axes', '')
                    if not isinstance(axes, str) or not axes:
                        nd = arr.ndim
                        if nd == 5:
                            axes = 'TCZYX'
                        elif nd == 4:
                            axes = 'CZYX'
                        elif nd == 3:
                            axes = 'ZYX'
                        elif nd == 2:
                            axes = 'YX'
                        else:
                            axes = 'C' * nd
                    axes = axes.upper()

                    if 'T' in axes:
                        t_axis = axes.index('T')
                        arr = np.take(arr, 0, axis=t_axis)
                        axes = axes.replace('T', '')

                    if 'A' in axes:
                        a_axis = axes.index('A')
                        arr = np.take(arr, 0, axis=a_axis)
                        axes = axes.replace('A', '')

                    if 'C' not in axes:
                        arr = arr[np.newaxis, ...]
                        axes = 'C' + axes
                    if 'Z' not in axes:
                        c_idx = axes.index('C')
                        arr = arr[:, np.newaxis, ...]
                        axes = axes[:c_idx+1] + 'Z' + axes[c_idx+1:]

                    perm = [axes.index(a) for a in 'CZYX']
                    arr = np.transpose(arr, perm)

                    return arr.astype(np.float32)
            except ImportError:
                self.ctx.log(self.ctx.t("log.warn_oiffile_missing"))
            except Exception as e:
                self.ctx.log(self.ctx.t("log.warn_oiffile_failed", err=e))

        if ext == '.oir':
            try:
                import oirfile
                with oirfile.OirFile(filepath) as oir:
                    arr = oir.asarray()
                    if arr.ndim == 2:
                        arr = arr[np.newaxis, np.newaxis, ...]
                    elif arr.ndim == 3:
                        arr = arr[np.newaxis, ...]
                    elif arr.ndim == 4:
                        if arr.shape[0] > 8 and arr.shape[1] <= 8:
                            arr = np.transpose(arr, (1, 0, 2, 3))
                    elif arr.ndim == 5:
                        arr = arr[0]
                        if arr.shape[0] > 8 and arr.shape[1] <= 8:
                            arr = np.transpose(arr, (1, 0, 2, 3))
                    return arr.astype(np.float32)
            except ImportError:
                self.ctx.log(self.ctx.t("log.warn_oirfile_missing"))
            except Exception as e:
                self.ctx.log(self.ctx.t("log.warn_oirfile_failed", err=e))

        if ext == '.nd2':
            try:
                import nd2
                with nd2.ND2File(filepath) as f:
                    arr = f.asarray()
                    if arr.ndim == 5:
                        arr = arr[0]
                    return arr.astype(np.float32)
            except ImportError:
                self.ctx.log(self.ctx.t("log.warn_nd2_missing"))
            except Exception as e:
                self.ctx.log(self.ctx.t("log.warn_nd2_failed", err=e))

        if ext == '.czi':
            try:
                import czifile
                with czifile.CziFile(filepath) as czi:
                    arr = czi.asarray()
                    axes = getattr(czi, 'axes', None)
                    if not isinstance(axes, str) or not axes:
                        nd = arr.ndim
                        if nd == 5:
                            axes = 'TCZYX'
                        elif nd == 4:
                            axes = 'CZYX'
                        elif nd == 3:
                            axes = 'ZYX'
                        elif nd == 2:
                            axes = 'YX'
                        else:
                            axes = 'C' * nd
                    axes = axes.upper()

                    if 'T' in axes:
                        t_axis = axes.index('T')
                        arr = np.take(arr, 0, axis=t_axis)
                        axes = axes.replace('T', '')

                    if 'S' in axes:
                        s_axis = axes.index('S')
                        arr = np.take(arr, 0, axis=s_axis)
                        axes = axes.replace('S', '')

                    if 'C' not in axes:
                        arr = arr[np.newaxis, ...]
                        axes = 'C' + axes
                    if 'Z' not in axes:
                        c_idx = axes.index('C')
                        arr = arr[:, np.newaxis, ...]
                        axes = axes[:c_idx+1] + 'Z' + axes[c_idx+1:]

                    perm = [axes.index(a) for a in 'CZYX']
                    arr = np.transpose(arr, perm)

                    return arr.astype(np.float32)
            except ImportError:
                self.ctx.log(self.ctx.t("log.warn_czifile_missing"))
                return None
            except Exception as e:
                self.ctx.log(self.ctx.t("log.warn_czifile_failed", err=e))
                return None

        if ext in ['.tif', '.tiff', '.ome.tif', '.ome.tiff']:
            try:
                with tifffile.TiffFile(filepath) as tif:
                    if len(tif.series) > 0:
                        series = tif.series[0]
                        arr = series.asarray()
                        axes = series.axes.upper()

                        if 'T' in axes and 'C' in axes:
                            t_axis = axes.index('T')
                            arr = np.take(arr, 0, axis=t_axis)
                            axes = axes.replace('T', '')

                        if axes == 'TZYX' and tif.imagej_metadata is not None:
                            ij = tif.imagej_metadata
                            n_channels = ij.get('channels', None)
                            n_frames = ij.get('frames', None)

                            if n_channels == 2:
                                axes = 'CZYX'
                            elif n_frames == 2:
                                axes = 'CZYX'
                            else:
                                arr = np.take(arr, 0, axis=0)
                                axes = 'ZYX'

                        if 'T' in axes:
                            t_axis = axes.index('T')
                            arr = np.take(arr, 0, axis=t_axis)
                            axes = axes.replace('T', '')

                        if axes == 'CZYX':
                            pass
                        elif axes == 'CYX':
                            arr = arr[:, np.newaxis, ...]
                        elif axes == 'ZYX':
                            arr = arr[np.newaxis, ...]
                        elif axes == 'ZCYX':
                            arr = np.transpose(arr, (1, 0, 2, 3))
                        elif axes == 'YX':
                            arr = arr[np.newaxis, np.newaxis, ...]
                        else:
                            self.ctx.log(self.ctx.t("log.warn_unknown_axes", axes=axes))
                            if arr.ndim == 4:
                                pass
                            elif arr.ndim == 3:
                                arr = arr[np.newaxis, ...]
                            elif arr.ndim == 2:
                                arr = arr[np.newaxis, np.newaxis, ...]

                        return arr.astype(np.float32)
                    else:
                        img = tif.asarray()
                        return img.astype(np.float32)
            except Exception as e:
                self.ctx.log(self.ctx.t("log.warn_tifffile_failed", err=e))
                try:
                    img = tifffile.imread(filepath)
                    return img.astype(np.float32)
                except Exception as e2:
                    self.ctx.log(self.ctx.t("log.warn_tifffile_imread_failed", err=e2))
                    return None

        self.ctx.log(self.ctx.t("log.warn_format_unsupported", ext=ext))
        return None

    def load_and_downsample(self, filepath):
        self.ctx.log(self.ctx.t("log.loading", name=os.path.basename(filepath)))

        try:
            img = self._read_image_universal(filepath)
            if img is None:
                return None

            if img.ndim == 2:
                img = img[np.newaxis, np.newaxis, ...]
            elif img.ndim == 3:
                if img.shape[0] <= 8 and img.shape[1] > img.shape[0]:
                    img = img[:, np.newaxis, ...]
                else:
                    img = img[np.newaxis, ...]
            elif img.ndim == 4:
                if img.shape[0] <= 8 and img.shape[1] > 4:
                    pass
                else:
                    img = np.transpose(img, (1, 0, 2, 3))
            elif img.ndim == 5:
                img = img[0]
                if img.shape[0] > 8 and img.shape[1] <= 8:
                    img = np.transpose(img, (1, 0, 2, 3))

            n_ch = img.shape[0]

            cell_idx = self.ctx.cell_channel_idx
            if cell_idx is None:
                cell_idx = 0
            if cell_idx < 0 or cell_idx >= n_ch:
                cell_idx = 0

            if self.ctx.has_dapi:
                dapi_idx = self.ctx.dapi_channel_idx
                if (dapi_idx is None or dapi_idx < 0 or
                        dapi_idx >= n_ch or dapi_idx == cell_idx):
                    alt = [k for k in range(n_ch) if k != cell_idx]
                    if alt:
                        dapi_idx = alt[0]
                    else:
                        self.ctx.has_dapi = False
                        dapi_idx = None

                if self.ctx.has_dapi:
                    img = img[[cell_idx, dapi_idx], ...]
                else:
                    img = img[cell_idx:cell_idx+1, ...]
            else:
                img = img[cell_idx:cell_idx+1, ...]

            if self.ctx.has_dapi and img.shape[0] < 2:
                self.ctx.log(self.ctx.t("log.warn_only_one_channel", n=img.shape[0]))
                img = np.concatenate([img, img], axis=0)

            c, z, y, x = img.shape

            self.ctx.is_2d = (z == 1)
            if self.ctx.is_2d:
                self.ctx.log(self.ctx.t("log.2d_image"))
            else:
                self.ctx.log(self.ctx.t("log.3d_stack", n=z))

            ds = max(1, int(self.ctx.downsample))
            if ds == 1:
                return img.astype(np.float32)

            ny, nx = y // ds, x // ds
            img_small = np.zeros((c, z, ny, nx), dtype=np.float32)

            for ch in range(c):
                for sl in range(z):
                    img_small[ch, sl] = resize(
                        img[ch, sl],
                        (ny, nx),
                        preserve_range=True,
                        anti_aliasing=True,
                        order=1
                    )

            return img_small

        except Exception as e:
            self.ctx.log(self.ctx.t("log.error", err=e))
            return None

    def get_file_paths(self):
        supported_exts = ['.tif', '.tiff', '.nd2', '.czi', '.lif', '.ome.tif', '.ome.tiff',
                          '.oib', '.oif', '.oir']
        if self.ctx.bridge is None:
            return []
        files = self.ctx.bridge._ask("file_paths", exts=supported_exts)
        return files or []

    def _extract_voxel_size_from_file(self, filepath):
        ext = os.path.splitext(filepath)[1].lower()
        result = {'xy_um': None, 'z_um': None, 'source': 'none'}

        if ext in ['.oib', '.oif']:
            try:
                import oiffile
                with oiffile.OifFile(filepath) as oib:
                    main = oib.mainfile
                    xy = None
                    z_um = None

                    try:
                        axis0 = main.get('Axis 0 Parameters', {}) or {}
                        px = axis0.get('Pixel Size')
                        if px is not None:
                            xy = float(px)
                    except Exception:
                        pass

                    try:
                        axis2 = main.get('Axis 2 Parameters', {}) or {}
                        pz = axis2.get('Pixel Size')
                        if pz is not None:
                            z_um = float(pz)
                    except Exception:
                        pass

                    if xy is not None and xy > 0:
                        result['xy_um'] = xy
                    if z_um is not None and z_um > 0:
                        result['z_um'] = z_um

                    if result['xy_um'] is not None or result['z_um'] is not None:
                        result['source'] = 'oib (oiffile)'
                        return result
            except ImportError:
                pass
            except Exception:
                pass

        if ext == '.oir':
            try:
                import oirfile
                with oirfile.OirFile(filepath) as oir:
                    xy = None
                    z_um = None

                    try:
                        scales = oir.coord_scales
                        if scales:
                            sx = scales.get('X')
                            sy = scales.get('Y')
                            sz = scales.get('Z')
                            if sx is not None and sy is not None:
                                xy = (float(sx) + float(sy)) / 2.0
                            elif sx is not None:
                                xy = float(sx)
                            if sz is not None:
                                z_um = float(sz)
                    except Exception:
                        pass

                    if xy is None or z_um is None:
                        try:
                            da = oir.asxarray()
                            if 'X' in da.coords and len(da.coords['X']) > 1:
                                xs = da.coords['X'].values
                                xy = float(abs(xs[1] - xs[0]))
                            if 'Z' in da.coords and len(da.coords['Z']) > 1:
                                zs = da.coords['Z'].values
                                z_um = float(abs(zs[1] - zs[0]))
                        except Exception:
                            pass

                    if xy is not None and xy > 0:
                        result['xy_um'] = xy
                    if z_um is not None and z_um > 0:
                        result['z_um'] = z_um

                    if result['xy_um'] is not None or result['z_um'] is not None:
                        result['source'] = 'oir (coord_scales)'
                        return result
            except ImportError:
                pass
            except Exception:
                pass

        if ext == '.nd2':
            try:
                import nd2
                with nd2.ND2File(filepath) as f:
                    vs = f.voxel_size()
                    if vs is not None:
                        result['xy_um'] = (vs.x + vs.y) / 2.0
                        result['z_um'] = vs.z
                        result['source'] = 'nd2 (voxel_size)'
                        return result
            except Exception:
                pass

        if ext == '.czi':
            try:
                import czifile
                with czifile.CziFile(filepath) as czi:
                    xml_str = None
                    try:
                        xml_str = czi.metadata()
                    except Exception:
                        pass

                    if xml_str is None:
                        return result

                    try:
                        root = ET.fromstring(xml_str)
                    except Exception:
                        return result

                    def _local(tag):
                        return tag.split('}')[-1] if isinstance(tag, str) else str(tag)

                    dist_xyz = {}

                    for el in root.iter():
                        if _local(el.tag) != 'Distance':
                            continue
                        dist_id = (el.attrib or {}).get('Id')
                        if dist_id not in ('X', 'Y', 'Z'):
                            continue
                        val = None
                        for child in el:
                            if _local(child.tag) == 'Value':
                                try:
                                    val = float(child.text)
                                except Exception:
                                    val = None
                                break
                        if val is not None:
                            dist_xyz[dist_id] = val

                    if 'X' in dist_xyz and 'Y' in dist_xyz:
                        xy_m = (abs(dist_xyz['X']) + abs(dist_xyz['Y'])) / 2.0
                        result['xy_um'] = xy_m * 1e6

                    if 'Z' in dist_xyz:
                        result['z_um'] = abs(dist_xyz['Z']) * 1e6

                    if result['xy_um'] is not None or result['z_um'] is not None:
                        result['source'] = 'czi (Distance)'
                        return result

                    for el in root.iter():
                        if _local(el.tag) != 'ImagePixelDistances':
                            continue
                        try:
                            parts = (el.text or '').split(',')
                            vals = [float(p.strip()) for p in parts if p.strip()]
                        except Exception:
                            vals = []
                        if len(vals) >= 2:
                            result['xy_um'] = (vals[0] + vals[1]) / 2.0
                            result['source'] = 'czi (ImagePixelDistances)'
                            break

                    if result['xy_um'] is not None:
                        return result

                    scaling_x = None
                    scaling_y = None
                    scaling_z = None

                    for el in root.iter():
                        tl = _local(el.tag)
                        if tl == 'ScalingX' and scaling_x is None:
                            try:
                                scaling_x = float(el.text)
                            except Exception:
                                pass
                        elif tl == 'ScalingY' and scaling_y is None:
                            try:
                                scaling_y = float(el.text)
                            except Exception:
                                pass
                        elif tl == 'ScalingZ' and scaling_z is None:
                            try:
                                scaling_z = float(el.text)
                            except Exception:
                                pass

                    if scaling_x is not None and scaling_y is not None:
                        if abs(scaling_x) < 1e-3:
                            xy = (abs(scaling_x) + abs(scaling_y)) / 2.0 * 1e6
                        else:
                            xy = (abs(scaling_x) + abs(scaling_y)) / 2.0
                        result['xy_um'] = xy

                    if scaling_z is not None:
                        if abs(scaling_z) < 1e-3:
                            z_um = abs(scaling_z) * 1e6
                        else:
                            z_um = abs(scaling_z)
                        result['z_um'] = z_um

                    if result['xy_um'] is not None or result['z_um'] is not None:
                        result['source'] = 'czi (Scaling)'
                        return result

                    return result
            except ImportError:
                pass
            except Exception:
                pass

        if ext in ['.tif', '.tiff', '.ome.tif', '.ome.tiff']:
            try:
                with tifffile.TiffFile(filepath) as tif:
                    if tif.ome_metadata is not None:
                        try:
                            root = ET.fromstring(tif.ome_metadata)
                            ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
                            pixels = root.find('.//ome:Pixels', ns)
                            if pixels is not None:
                                px = pixels.get('PhysicalSizeX')
                                py = pixels.get('PhysicalSizeY')
                                pz = pixels.get('PhysicalSizeZ')
                                unit_x = pixels.get('PhysicalSizeXUnit', 'µm')
                                unit_y = pixels.get('PhysicalSizeYUnit', 'µm')
                                unit_z = pixels.get('PhysicalSizeZUnit', 'µm')

                                def to_um(val, unit):
                                    if val is None:
                                        return None
                                    try:
                                        v = float(val)
                                    except Exception:
                                        return None
                                    unit = (unit or 'µm').lower().strip()
                                    if 'nm' in unit:
                                        return v / 1000.0
                                    if unit in ('µm', 'um', 'micron', 'microns',
                                                'micrometer', 'micrometers'):
                                        return v
                                    if 'mm' in unit:
                                        return v * 1000.0
                                    if 'cm' in unit:
                                        return v * 10000.0
                                    if unit == 'm':
                                        return v * 1e6
                                    return v

                                xy = None
                                if px is not None and py is not None:
                                    xy = (to_um(px, unit_x) + to_um(py, unit_y)) / 2.0
                                result['xy_um'] = xy
                                result['z_um'] = to_um(pz, unit_z) if pz is not None else None
                                result['source'] = 'ome-xml'
                                return result
                        except Exception:
                            pass

                    if tif.imagej_metadata is not None:
                        ij = tif.imagej_metadata

                        z_spacing = ij.get('spacing', None)
                        unit_ij = ij.get('unit', 'µm')

                        z_um = None
                        if z_spacing is not None:
                            try:
                                z_um = float(z_spacing)
                                u = (unit_ij or '').lower()
                                if 'nm' in u:
                                    z_um = z_um / 1000.0
                                elif 'mm' in u:
                                    z_um = z_um * 1000.0
                            except Exception:
                                z_um = None

                        if z_um is None:
                            info_str = ij.get('Info', '') or ''
                            m = re.search(r'spacing\s*=\s*([0-9.eE+\-]+)', info_str)
                            if m:
                                try:
                                    z_um = float(m.group(1))
                                except Exception:
                                    pass

                        if z_um is None:
                            info_str = ij.get('Info', '') or ''
                            m = re.search(r'z[_\- ]?step\s*=\s*([0-9.eE+\-]+)',
                                          info_str, re.IGNORECASE)
                            if m:
                                try:
                                    z_um = float(m.group(1))
                                except Exception:
                                    pass

                        xy = None
                        try:
                            tags = tif.pages[0].tags
                            xr = tags.get('XResolution')
                            yr = tags.get('YResolution')
                            res_unit_tag = tags.get('ResolutionUnit')

                            res_unit = 2
                            if res_unit_tag is not None:
                                res_unit = res_unit_tag.value

                            xr_num, xr_den = (xr.value if xr else (0, 1))
                            yr_num, yr_den = (yr.value if yr else (0, 1))
                            x_ppu = xr_num / xr_den if xr_den != 0 else 0
                            y_ppu = yr_num / yr_den if yr_den != 0 else 0

                            if x_ppu > 0 and y_ppu > 0:
                                if res_unit == 2:
                                    x_um = 25400.0 / x_ppu
                                    y_um = 25400.0 / y_ppu
                                elif res_unit == 3:
                                    x_um = 10000.0 / x_ppu
                                    y_um = 10000.0 / y_ppu
                                else:
                                    x_um = 1.0 / x_ppu
                                    y_um = 1.0 / y_ppu

                                xy = (x_um + y_um) / 2.0
                        except Exception:
                            pass

                        result['xy_um'] = xy
                        result['z_um'] = z_um
                        result['source'] = 'imagej'
                        return result

                    try:
                        tags = tif.pages[0].tags
                        xr = tags.get('XResolution')
                        yr = tags.get('YResolution')
                        res_unit_tag = tags.get('ResolutionUnit')

                        res_unit = res_unit_tag.value if res_unit_tag is not None else 2

                        if xr is not None and yr is not None:
                            xr_num, xr_den = xr.value
                            yr_num, yr_den = yr.value
                            x_ppu = xr_num / xr_den if xr_den != 0 else 0
                            y_ppu = yr_num / yr_den if yr_den != 0 else 0
                            if x_ppu > 0 and y_ppu > 0:
                                if res_unit == 2:
                                    x_um = 25400.0 / x_ppu
                                    y_um = 25400.0 / y_ppu
                                elif res_unit == 3:
                                    x_um = 10000.0 / x_ppu
                                    y_um = 10000.0 / y_ppu
                                else:
                                    x_um = 1.0 / x_ppu
                                    y_um = 1.0 / y_ppu
                                xy = (x_um + y_um) / 2.0
                                result['xy_um'] = xy
                                result['source'] = 'tiff-tags'
                                return result
                    except Exception:
                        pass

            except Exception:
                pass

        return result

    def extract_and_confirm_metadata(self, filepath):
        self.ctx.log(self.ctx.t("log.voxel_size_title"))
        self.ctx.log(self.ctx.t("log.file_label", name=os.path.basename(filepath)))

        meta = self._extract_voxel_size_from_file(filepath)
        self.ctx._gui_metadata_cache = meta

        if getattr(self.ctx, "_gui_scale_mode", "auto") == "manual":
            self.ctx.log(self.ctx.t("log.mode_manual"))
            self.ctx.log(self.ctx.t("log.scale_value",
                                    xy=self.ctx.xy_scale_um,
                                    z=self.ctx.z_scale_um))
            self.ctx.scale_initialized = True
            self.ctx._input_xy_um = self.ctx.xy_scale_um
            self.ctx._input_z_um = self.ctx.z_scale_um
            return

        if meta['source'] == 'none' or (meta['xy_um'] is None and meta['z_um'] is None):
            self.ctx.log(self.ctx.t("log.no_voxel_metadata"))
            self.ctx.log(self.ctx.t("log.use_panel_settings"))
            if self.ctx.xy_scale_um is None:
                self.ctx.xy_scale_um = 0.144
            if self.ctx.z_scale_um is None:
                self.ctx.z_scale_um = 1.0
            self.ctx.scale_initialized = True
            self.ctx._input_xy_um = self.ctx.xy_scale_um
            self.ctx._input_z_um = self.ctx.z_scale_um
            self.ctx.log(self.ctx.t("log.scale_value",
                                    xy=self.ctx.xy_scale_um,
                                    z=self.ctx.z_scale_um))
            return

        self.ctx.log(self.ctx.t("log.source_label", src=meta['source']))
        if meta['xy_um'] is not None:
            self.ctx.log(self.ctx.t("log.xy_value", val=meta['xy_um']))
        else:
            self.ctx.log(self.ctx.t("log.xy_not_defined"))
        if meta['z_um'] is not None:
            self.ctx.log(self.ctx.t("log.z_value", val=meta['z_um']))
        else:
            self.ctx.log(self.ctx.t("log.z_not_defined"))

        if self.ctx.xy_scale_um is None:
            self.ctx.xy_scale_um = 0.144
        if self.ctx.z_scale_um is None:
            self.ctx.z_scale_um = 1.0

        if meta['xy_um'] is not None and meta['xy_um'] > 0:
            self.ctx.xy_scale_um = meta['xy_um']
        if meta['z_um'] is not None and meta['z_um'] > 0:
            self.ctx.z_scale_um = meta['z_um']

        self.ctx.scale_initialized = True
        self.ctx._input_xy_um = self.ctx.xy_scale_um
        self.ctx._input_z_um = self.ctx.z_scale_um
        self.ctx.log(self.ctx.t("log.scale_value",
                                xy=self.ctx.xy_scale_um,
                                z=self.ctx.z_scale_um))

    def get_scale_from_user(self):
        pass

class Optics:
    def __init__(self, ctx, file_io):
        self.ctx = ctx
        self.file_io = file_io

    def _extract_optics_from_file(self, filepath):
        ext = os.path.splitext(filepath)[1].lower()
        result = {
            'channels': [],
            'dapi_channel_index': None,
            'astro_channel_index': None,
            'source': 'none'
        }

        if ext == '.nd2':
            try:
                import nd2
                with nd2.ND2File(filepath) as f:
                    md = f.metadata
                    channels = getattr(md, 'channels', None) or []

                    for i, ch_wrapper in enumerate(channels[:2]):
                        ch_meta = getattr(ch_wrapper, 'channel', None)
                        microscope = getattr(ch_wrapper, 'microscope', None)

                        ch_info = {
                            'index': i,
                            'name': None,
                            'emission_um': None,
                            'excitation_um': None,
                            'na': None,
                            'magnification': None
                        }

                        if ch_meta is not None:
                            ch_info['name'] = getattr(ch_meta, 'name', None)

                            em_nm = getattr(ch_meta, 'emissionLambdaNm', None)
                            if em_nm is not None:
                                ch_info['emission_um'] = float(em_nm) / 1000.0

                            ex_nm = getattr(ch_meta, 'excitationLambdaNm', None)
                            if ex_nm is not None:
                                ch_info['excitation_um'] = float(ex_nm) / 1000.0

                        if microscope is not None:
                            na = getattr(microscope, 'objectiveNumericalAperture', None)
                            if na is not None:
                                ch_info['na'] = float(na)
                            mag = getattr(microscope, 'objectiveMagnification', None)
                            if mag is not None:
                                ch_info['magnification'] = float(mag)

                        result['channels'].append(ch_info)

                    self._assign_dapi_by_metadata(result)

                    if result['channels']:
                        result['source'] = 'nd2'
                    return result
            except Exception:
                pass

        if ext == '.czi':
            try:
                import czifile
                with czifile.CziFile(filepath) as czi:
                    try:
                        xml_str = czi.metadata()
                    except Exception:
                        xml_str = None

                    if xml_str is None:
                        return result

                    try:
                        root = ET.fromstring(xml_str)
                    except Exception:
                        return result

                    def _local(tag):
                        return tag.split('}')[-1] if isinstance(tag, str) else str(tag)

                    na_global = None
                    mag_global = None

                    for el in root.iter():
                        if _local(el.tag) != 'Objective':
                            continue
                        for attr_name in ('LensNA', 'NumericalAperture'):
                            attr_v = el.get(attr_name)
                            if attr_v is not None:
                                try:
                                    na_global = float(attr_v)
                                    break
                                except Exception:
                                    pass
                        for attr_name in ('NominalMagnification', 'Magnification'):
                            attr_v = el.get(attr_name)
                            if attr_v is not None:
                                try:
                                    mag_global = float(attr_v)
                                    break
                                except Exception:
                                    pass

                    channels_found = []
                    for idx, el in enumerate(root.iter()):
                        if _local(el.tag) != 'Channel':
                            continue
                        if idx >= 2:
                            break
                        ch_info = {
                            'index': len(channels_found),
                            'name': el.get('Name') or el.get('Id'),
                            'emission_um': None,
                            'excitation_um': None,
                            'na': na_global,
                            'magnification': mag_global
                        }
                        for attr_name in ('EmissionWavelength', 'EmissionWavelengthNm'):
                            v = el.get(attr_name)
                            if v is not None:
                                try:
                                    val = float(v)
                                    if val > 10:
                                        val = val / 1000.0
                                    ch_info['emission_um'] = val
                                    break
                                except Exception:
                                    pass
                        for attr_name in ('ExcitationWavelength', 'ExcitationWavelengthNm'):
                            v = el.get(attr_name)
                            if v is not None:
                                try:
                                    val = float(v)
                                    if val > 10:
                                        val = val / 1000.0
                                    ch_info['excitation_um'] = val
                                    break
                                except Exception:
                                    pass
                        channels_found.append(ch_info)
                        if len(channels_found) >= 2:
                            break

                    result['channels'] = channels_found
                    self._assign_dapi_by_metadata(result)

                    if result['channels']:
                        result['source'] = 'czi'
                    return result
            except ImportError:
                pass
            except Exception:
                pass

        if ext in ['.tif', '.tiff', '.ome.tif', '.ome.tiff']:
            try:
                with tifffile.TiffFile(filepath) as tif:
                    if tif.ome_metadata is not None:
                        root = ET.fromstring(tif.ome_metadata)
                        ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

                        objective = root.find('.//ome:Objective', ns)
                        na_global = None
                        mag_global = None
                        if objective is not None:
                            na_attr = objective.get('LensNA')
                            if na_attr is not None:
                                na_global = float(na_attr)
                            mag_attr = objective.get('NominalMagnification')
                            if mag_attr is not None:
                                mag_global = float(mag_attr)

                        channels = root.findall('.//ome:Channel', ns)
                        for i, ch_el in enumerate(channels[:2]):
                            ch_info = {
                                'index': i,
                                'name': ch_el.get('Name') or ch_el.get('ID'),
                                'emission_um': None,
                                'excitation_um': None,
                                'na': na_global,
                                'magnification': mag_global
                            }

                            em = ch_el.get('EmissionWavelength')
                            em_unit = ch_el.get('EmissionWavelengthUnit', 'nm')
                            if em is not None:
                                em = float(em)
                                if 'nm' in em_unit.lower():
                                    em = em / 1000.0
                                ch_info['emission_um'] = em

                            ex = ch_el.get('ExcitationWavelength')
                            ex_unit = ch_el.get('ExcitationWavelengthUnit', 'nm')
                            if ex is not None:
                                ex = float(ex)
                                if 'nm' in ex_unit.lower():
                                    ex = ex / 1000.0
                                ch_info['excitation_um'] = ex

                            result['channels'].append(ch_info)

                        self._assign_dapi_by_metadata(result)

                        if result['channels']:
                            result['source'] = 'ome-xml'
                        return result
            except Exception:
                pass

        return result

    def _assign_dapi_by_metadata(self, result):
        channels = result['channels']
        if len(channels) < 2:
            return

        dapi_name_keywords = ['dapi', 'dna', 'hoechst', 'blue', ' nucleus', 'nuclei']
        for ch in channels:
            name = (ch.get('name') or '').lower()
            for kw in dapi_name_keywords:
                if kw in name:
                    result['dapi_channel_index'] = ch['index']
                    result['astro_channel_index'] = 1 - ch['index']
                    return

        for ch in channels:
            ex = ch.get('excitation_um')
            if ex is not None and 0.350 <= ex <= 0.410:
                result['dapi_channel_index'] = ch['index']
                result['astro_channel_index'] = 1 - ch['index']
                return

        for ch in channels:
            em = ch.get('emission_um')
            if em is not None and 0.440 <= em <= 0.500:
                result['dapi_channel_index'] = ch['index']
                result['astro_channel_index'] = 1 - ch['index']
                return

    def determine_channels(self, filepath):
        self.ctx.log(self.ctx.t("log.channels_title"))

        optics = self._extract_optics_from_file(filepath)
        channels = optics.get('channels', [])

        if channels:
            self.ctx.log(self.ctx.t("log.source_label", src=optics['source']))
            for ch in channels:
                if ch.get('name'):
                    self.ctx.log(self.ctx.t("log.channel_line",
                                            idx=ch['index'],
                                            name=ch['name'],
                                            em=(ch.get('emission_um') or 0) * 1000,
                                            ex=(ch.get('excitation_um') or 0) * 1000))
                else:
                    self.ctx.log(self.ctx.t("log.channel_line_no_name",
                                            idx=ch['index'],
                                            em=(ch.get('emission_um') or 0) * 1000,
                                            ex=(ch.get('excitation_um') or 0) * 1000))

        auto_cell_idx = optics.get('astro_channel_index')
        auto_dapi_idx = optics.get('dapi_channel_index')

        img = self.file_io._read_image_universal(filepath)
        if img is None:
            self.ctx.log(self.ctx.t("log.file_read_failed"))
            self.ctx.has_dapi = False
            self.ctx.cell_channel_idx = 0
            self.ctx.dapi_channel_idx = None
            self.ctx.selected_channels = [0]
            self.ctx.astrocyte_channel = 0
            self.ctx.dapi_channel = 0
            return

        if img.ndim == 2:
            img = img[np.newaxis, np.newaxis, ...]
        elif img.ndim == 3:
            if img.shape[0] <= 8 and img.shape[1] > img.shape[0]:
                img = img[:, np.newaxis, ...]
            else:
                img = img[np.newaxis, ...]
        elif img.ndim == 4:
            if img.shape[0] <= 8:
                pass
            else:
                img = np.transpose(img, (1, 0, 2, 3))
        elif img.ndim == 5:
            img = img[0]
            if img.shape[0] > 8 and img.shape[1] <= 8:
                img = np.transpose(img, (1, 0, 2, 3))

        n_ch = img.shape[0]
        self.ctx.n_channels_original = int(n_ch)
        self.ctx.log(self.ctx.t("log.channels_found", n=n_ch))

        if n_ch == 1:
            self.ctx.has_dapi = False
            self.ctx.cell_channel_idx = 0
            self.ctx.dapi_channel_idx = None
            self.ctx.selected_channels = [0]
            self.ctx.astrocyte_channel = 0
            self.ctx.dapi_channel = 0
            self.ctx.log(self.ctx.t("log.single_channel"))
            return

        cell_default = 0
        if auto_cell_idx is not None and 0 <= auto_cell_idx < n_ch:
            cell_default = int(auto_cell_idx)

        dapi_default = None
        if auto_dapi_idx is not None and 0 <= auto_dapi_idx < n_ch:
            dapi_default = int(auto_dapi_idx)

        if dapi_default == cell_default:
            dapi_default = None

        has_dapi_default = (dapi_default is not None)

        def norm(m):
            v = m[m > 0]
            if len(v) == 0:
                return m
            p0, p99 = np.percentile(v, [1, 99])
            return np.clip((m - p0) / (p99 - p0 + 1e-6), 0, 1)

        mips = [norm(np.max(img[i], axis=0)) for i in range(n_ch)]

        self.ctx.log(self.ctx.t("log.asking_channel_choice"))
        res = self.ctx.bridge._ask(
            "channel_selector",
            mips=mips,
            cell_idx_default=cell_default,
            dapi_idx_default=dapi_default,
            has_dapi_default=has_dapi_default,
        )

        if res is None:
            self.ctx.log(self.ctx.t("log.channel_cancelled"))
            self.ctx.cell_channel_idx = cell_default
            self.ctx.has_dapi = has_dapi_default
            self.ctx.dapi_channel_idx = dapi_default
            if has_dapi_default:
                self.ctx.selected_channels = [cell_default, dapi_default]
            else:
                self.ctx.selected_channels = [cell_default]
        else:
            self.ctx.selected_channels = list(res["selected_channels"])
            self.ctx.cell_channel_idx = int(res["cell_idx"])
            self.ctx.has_dapi = bool(res["has_dapi"])
            self.ctx.dapi_channel_idx = (int(res["dapi_idx"])
                                         if res["dapi_idx"] is not None else None)

        self.ctx.astrocyte_channel = 0
        self.ctx.dapi_channel = 1 if self.ctx.has_dapi else 0

        self.ctx.log(self.ctx.t("log.cells_channel", idx=self.ctx.cell_channel_idx))
        if self.ctx.has_dapi:
            self.ctx.log(self.ctx.t("log.dapi_channel", idx=self.ctx.dapi_channel_idx))
        else:
            self.ctx.log(self.ctx.t("log.dapi_not_used"))
        self.ctx.log(self.ctx.t("log.channels_in_analysis",
                                ch=self.ctx.selected_channels))

    def extract_and_confirm_optics(self, filepath):
        self.ctx.log(self.ctx.t("log.optics_title"))
        self.ctx.log(self.ctx.t("log.file_label", name=os.path.basename(filepath)))

        optics = self._extract_optics_from_file(filepath)
        channels = optics.get('channels', [])

        if getattr(self.ctx, "_gui_deconv_mode", "auto") == "manual":
            self.ctx.log(self.ctx.t("log.optics_manual"))
            self.ctx.log(self.ctx.t("log.optics_cells",
                                    na=self.ctx._deconv_params['na'],
                                    wl=self.ctx._deconv_params['emission_wavelength_um'] * 1000))
            self.ctx.log(self.ctx.t("log.optics_dapi",
                                    na=self.ctx._dapi_deconv_params['na'],
                                    wl=self.ctx._dapi_deconv_params['emission_wavelength_um'] * 1000))
            return

        if not channels:
            self.ctx.log(self.ctx.t("log.optics_no_metadata"))
            self.ctx.log(self.ctx.t("log.use_panel_settings"))
            self.ctx.log(self.ctx.t("log.optics_cells",
                                    na=self.ctx._deconv_params['na'],
                                    wl=self.ctx._deconv_params['emission_wavelength_um'] * 1000))
            self.ctx.log(self.ctx.t("log.optics_dapi",
                                    na=self.ctx._dapi_deconv_params['na'],
                                    wl=self.ctx._dapi_deconv_params['emission_wavelength_um'] * 1000))
            return

        self.ctx.log(self.ctx.t("log.source_label", src=optics['source']))

        dapi_idx = self.ctx.dapi_channel
        astro_idx = self.ctx.astrocyte_channel

        astro_ch = channels[astro_idx] if astro_idx < len(channels) else {}
        dapi_ch = channels[dapi_idx] if dapi_idx < len(channels) else {}

        if astro_ch.get('na') is not None and astro_ch['na'] > 0:
            self.ctx._deconv_params['na'] = astro_ch['na']
        if astro_ch.get('emission_um') is not None and astro_ch['emission_um'] > 0:
            self.ctx._deconv_params['emission_wavelength_um'] = astro_ch['emission_um']

        if dapi_ch.get('na') is not None and dapi_ch['na'] > 0:
            self.ctx._dapi_deconv_params['na'] = dapi_ch['na']
        else:
            self.ctx._dapi_deconv_params['na'] = self.ctx._deconv_params['na']
        if dapi_ch.get('emission_um') is not None and dapi_ch['emission_um'] > 0:
            self.ctx._dapi_deconv_params['emission_wavelength_um'] = dapi_ch['emission_um']

        self.ctx.log(self.ctx.t("log.optics_cells",
                                na=self.ctx._deconv_params['na'],
                                wl=self.ctx._deconv_params['emission_wavelength_um'] * 1000))
        self.ctx.log(self.ctx.t("log.optics_dapi",
                                na=self.ctx._dapi_deconv_params['na'],
                                wl=self.ctx._dapi_deconv_params['emission_wavelength_um'] * 1000))

class Deconvolution:
    def __init__(self, ctx):
        self.ctx = ctx

    def _deconvolve_3d_stack(self, image_stack, pixel_xy_um, pixel_z_um):
        self.ctx.log(self.ctx.t("log.deconv_start"))

        z_dim, y_dim, x_dim = image_stack.shape

        sigma_xy_bg = 2.0
        sigma_z_bg = 5.0
        radius_xy_bg = int(sigma_xy_bg * 3) + 1
        radius_z_bg = int(sigma_z_bg * 3) + 1

        z = np.arange(-radius_z_bg, radius_z_bg + 1)
        y = np.arange(-radius_xy_bg, radius_xy_bg + 1)
        x = np.arange(-radius_xy_bg, radius_xy_bg + 1)
        Z, Y, X = np.meshgrid(z, y, x, indexing='ij')

        psf_kernel_bg = np.exp(-((X**2 + Y**2) / (2 * sigma_xy_bg**2) + (Z**2) / (2 * sigma_z_bg**2)))
        psf_kernel_bg /= np.sum(psf_kernel_bg)

        img_bg = convolve(image_stack, psf_kernel_bg, mode='reflect')

        img_subtracted = image_stack - img_bg * 1.0
        img_subtracted = np.clip(img_subtracted, 0, None)

        sigma_unsharp = 0.3
        img_blur_small = gaussian_filter(img_subtracted, sigma=sigma_unsharp)

        unsharp_weight = 0.4
        img_sharp = img_subtracted + unsharp_weight * (img_subtracted - img_blur_small)
        img_sharp = np.clip(img_sharp, 0, None)

        img_deconv = median_filter(
            img_sharp,
            size=(1, 2, 2),
            mode='reflect'
        )

        original_sum = np.sum(image_stack)
        if original_sum > 0 and np.sum(img_deconv) > 0:
            scale_factor = original_sum / (np.sum(img_deconv) + 1e-6)
            img_deconv = img_deconv * scale_factor

        to_delete = self.ctx.bridge._ask(
            "slice_selector",
            stack=img_deconv,
            title=self.ctx.t("dlg.slices.title_selector"),
        )

        if to_delete and len(to_delete) > 0:
            self.ctx.log(self.ctx.t("log.deconv_slices_deleted",
                                    n=len(to_delete), z=to_delete))
            img_deconv[to_delete, :, :] = 0.0
        else:
            self.ctx.log(self.ctx.t("log.deconv_all_slices_kept"))

        self.ctx.log(self.ctx.t("log.deconv_done", shape=img_deconv.shape))
        return img_deconv.astype(np.float32)

    def _deconvolve_3d_stack_dapi(self, image_stack, pixel_xy_um, pixel_z_um):
        self.ctx.log(self.ctx.t("log.deconv_dapi_start"))

        z_dim, y_dim, x_dim = image_stack.shape

        na = self.ctx._dapi_deconv_params.get('na', 0.95)
        emission_um = self.ctx._dapi_deconv_params.get('emission_wavelength_um', 0.460)

        if na > 0 and emission_um > 0:
            airy_radius_um = 0.61 * emission_um / na
            sigma_xy_um = airy_radius_um / 2.5
            sigma_z_um = airy_radius_um / 1.5

            sigma_xy = sigma_xy_um / pixel_xy_um if pixel_xy_um > 0 else 0.6
            sigma_z = sigma_z_um / pixel_z_um if pixel_z_um > 0 else 1.5
        else:
            sigma_xy = 0.6
            sigma_z = 1.5

        sigma_xy = max(0.3, min(sigma_xy, 3.0))
        sigma_z = max(0.5, min(sigma_z, 5.0))

        img_blur = gaussian_filter(image_stack, sigma=(sigma_z, sigma_xy, sigma_xy), mode='reflect')

        weight = 1.2
        img_sharp = image_stack + weight * (image_stack - img_blur)
        img_sharp = np.clip(img_sharp, 0, None)

        original_sum = np.sum(image_stack)
        if original_sum > 0:
            scale_factor = original_sum / (np.sum(img_sharp) + 1e-6)
            img_deconv = img_sharp * scale_factor
        else:
            img_deconv = img_sharp

        img_deconv = median_filter(img_deconv, size=(3, 3, 3), mode='reflect')

        to_delete = self.ctx.bridge._ask(
            "slice_selector",
            stack=img_deconv,
            title=self.ctx.t("dlg.slices.title_selector"),
        )

        if to_delete and len(to_delete) > 0:
            self.ctx.log(self.ctx.t("log.deconv_slices_deleted",
                                    n=len(to_delete), z=to_delete))
            img_deconv[to_delete, :, :] = 0.0
        else:
            self.ctx.log(self.ctx.t("log.deconv_all_slices_kept"))

        self.ctx.log(self.ctx.t("log.deconv_dapi_done", shape=img_deconv.shape))
        return img_deconv.astype(np.float32)

class Segmentation:
    def __init__(self, ctx):
        self.ctx = ctx

    def measure_background(self, astrocyte_stack):
        finite_mask = np.isfinite(astrocyte_stack)
        if not np.any(finite_mask):
            return 0.0

        clean = astrocyte_stack[finite_mask]

        dark_percentile = np.percentile(clean, self.ctx.background_percentile)
        dark_mask = clean <= dark_percentile

        if not np.any(dark_mask):
            background = float(np.min(clean))
        else:
            background = float(np.median(clean[dark_mask]))

        if not np.isfinite(background):
            background = 0.0

        return background

    def manual_soma_setting(self, astrocyte_stack):
        self.ctx.log(self.ctx.t("log.soma_title"))

        z_dim, y_dim, x_dim = astrocyte_stack.shape
        astro = astrocyte_stack.astype(np.float32)
        astro_smooth = gaussian_filter(astro, sigma=0.5)
        p0, p99 = np.percentile(astro_smooth[astro_smooth > 0], [1, 99])
        astro_norm = np.clip((astro_smooth - p0) / (p99 - p0 + 1e-6), 0, 1)
        mip = np.max(astro_norm, axis=0)

        points = self.ctx.bridge._ask("soma_picker", mip=mip,
                                      astro_norm=astro_norm)

        if not points:
            self.ctx.log(self.ctx.t("log.soma_no_points"))
            return np.zeros_like(astrocyte_stack, dtype=np.int32)

        self.ctx.log(self.ctx.t("log.soma_points_set", n=len(points)))

        bfs_seeds = np.zeros_like(astrocyte_stack, dtype=np.int32)
        click_centers = []

        for i, pt in enumerate(points, 1):
            if len(pt) == 3:
                py, px, pz = pt
            else:
                py, px = pt
                pz = None

            soma_py = int(round(py))
            soma_px = int(round(px))
            soma_py = max(0, min(y_dim - 1, soma_py))
            soma_px = max(0, min(x_dim - 1, soma_px))

            if pz is not None:
                best_z_click = int(pz)
                best_z_click = max(0, min(z_dim - 1, best_z_click))
            else:
                best_z_click = 0
                best_score_click = -float('inf')
                for z in range(z_dim):
                    window = astro_norm[z,
                    max(0, soma_py - 2):min(y_dim, soma_py + 3),
                    max(0, soma_px - 2):min(x_dim, soma_px + 3)]
                    if window.size > 0:
                        score = np.mean(window)
                        if score > best_score_click:
                            best_score_click = score
                            best_z_click = z

            bfs_seeds[best_z_click, soma_py, soma_px] = i
            click_centers.append((best_z_click, soma_py, soma_px))

        ok = self.ctx.bridge._ask(
            "soma_confirm",
            astro_norm=astro_norm,
            click_centers=click_centers,
        )
        if not ok:
            self.ctx.log(self.ctx.t("log.soma_cancelled"))
            return np.zeros_like(astrocyte_stack, dtype=np.int32)

        self.ctx.log(self.ctx.t("log.soma_starting"))

        neighbors_26 = [(dz, dy, dx) for dz in (-1,0,1) for dy in (-1,0,1) for dx in (-1,0,1) if not (dz==0 and dy==0 and dx==0)]
        soma_segmented_3d = np.zeros_like(astrocyte_stack, dtype=np.int32)

        for i, (cz, cy, cx) in enumerate(click_centers, 1):
            if not (0 <= cz < z_dim and 0 <= cy < y_dim and 0 <= cx < x_dim):
                continue

            queue = deque()
            soma_segmented_3d[cz, cy, cx] = i
            queue.append((cz, cy, cx))

            initial_intensity = astro_norm[cz, cy, cx]
            intensity_threshold = max(0.15, initial_intensity * 0.70)
            max_radius_xy = 12
            if self.ctx.scale_initialized:
                xy_um = self.ctx.xy_scale_um * self.ctx.downsample
                max_radius_xy = min(12, int(7.0 / xy_um))
            max_radius_z = max(2, int(max_radius_xy * 0.5))

            while queue:
                z, y, x = queue.popleft()

                if ((z - cz) / max_radius_z)**2 + ((y - cy) / max_radius_xy)**2 + ((x - cx) / max_radius_xy)**2 > 1.0:
                    continue

                for dz, dy, dx in neighbors_26:
                    nz, ny, nx = z + dz, y + dy, x + dx
                    if not (0 <= nz < z_dim and 0 <= ny < y_dim and 0 <= nx < x_dim):
                        continue
                    if soma_segmented_3d[nz, ny, nx] != 0:
                        continue

                    if astro_norm[nz, ny, nx] >= intensity_threshold:
                        soma_segmented_3d[nz, ny, nx] = i
                        queue.append((nz, ny, nx))

        self.ctx._soma_labels = soma_segmented_3d

        self.ctx.log(self.ctx.t("log.soma_done"))
        return soma_segmented_3d.astype(np.int32)

    def segment_nuclei_local_from_soma(self, soma_mask, dapi_stack):
        self.ctx.log(self.ctx.t("log.nuclei_title"))

        if not getattr(self.ctx, "has_dapi", True):
            self.ctx.log(self.ctx.t("log.nuclei_no_dapi"))
            return np.zeros_like(dapi_stack, dtype=np.int32)

        clean_nuclei = np.zeros_like(dapi_stack, dtype=np.int32)
        unique_ids = np.unique(soma_mask)[1:]
        if len(unique_ids) == 0:
            self.ctx.log(self.ctx.t("log.nuclei_no_soma"))
            return clean_nuclei

        for cid in unique_ids:
            soma_mask_cid = (soma_mask == cid)
            dilated_soma = binary_dilation(soma_mask_cid, iterations=3)

            dapi_crop = dapi_stack[dilated_soma]
            if len(dapi_crop) < 10: continue

            dapi_vals = dapi_crop[dapi_crop > 0]
            if len(dapi_vals) < 10: continue

            try:
                thresh = threshold_otsu(dapi_vals)
                if np.isnan(thresh) or thresh == 0:
                    thresh = np.percentile(dapi_vals, 40)
            except:
                thresh = np.percentile(dapi_vals, 40)

            nuc_candidate_mask = np.zeros_like(dapi_stack, dtype=bool)
            nuc_candidate_mask[dilated_soma] = dapi_crop >= thresh

            selem = ball(1)
            nuc_candidate_mask = binary_opening(nuc_candidate_mask, structure=selem)
            nuc_candidate_mask = remove_small_objects(nuc_candidate_mask, min_size=15)

            if np.sum(nuc_candidate_mask) == 0: continue

            distance = distance_transform_edt(nuc_candidate_mask)
            peak = np.unravel_index(np.argmax(distance), distance.shape)
            markers = np.zeros_like(distance, dtype=int)
            markers[peak] = 1

            try:
                watershed_result = watershed(-distance, markers, mask=nuc_candidate_mask)
                final_nuc_mask = (watershed_result == 1)
            except:
                labeled_nuc, num_nuc = label(nuc_candidate_mask.astype(int),
                                             structure=np.ones((3,3,3)))
                if num_nuc > 0:
                    sizes = [np.sum(labeled_nuc == i) for i in range(1, num_nuc+1)]
                    largest_label = np.argmax(sizes) + 1
                    final_nuc_mask = (labeled_nuc == largest_label)
                else:
                    final_nuc_mask = nuc_candidate_mask

            clean_nuclei[final_nuc_mask] = cid

        self.ctx.log(self.ctx.t("log.nuclei_found",
                                n=len(np.unique(clean_nuclei)) - 1))
        return clean_nuclei

    def _grow_cells_by_components(self, astro_norm, soma_mask, growth_mask,
                                  unique_ids, dip_hard_threshold=0.5):
        import heapq
        from scipy.ndimage import laplace

        z_dim, y_dim, x_dim = astro_norm.shape

        astro = astro_norm.astype(np.float32)

        sigma_xy_bg = max(30, min(y_dim, x_dim) // 15)
        local_mean = uniform_filter(
            astro,
            size=(3, sigma_xy_bg, sigma_xy_bg),
            mode='reflect'
        )
        local_mean = np.maximum(local_mean, 1e-6)

        dip_mask = astro < (local_mean * 0.5)

        soma_bool = (soma_mask > 0)

        graph_mask = (growth_mask | soma_bool) & (~dip_mask)
        graph_mask = graph_mask | soma_bool

        if not np.any(graph_mask):
            self.ctx.log(self.ctx.t("log.cells_no_growth"))
            return soma_mask.astype(np.int32)

        labeled, num_comp = label(graph_mask, structure=np.ones((3, 3, 3)))

        component_somas = {}
        soma_coords_by_cid = {}

        for cid in unique_ids:
            soma_mask_cid = (soma_mask == cid)
            if not np.any(soma_mask_cid):
                continue
            coords = np.argwhere(soma_mask_cid)
            soma_coords_by_cid[cid] = coords

            comp_ids = np.unique(labeled[soma_mask_cid])
            comp_ids = comp_ids[comp_ids > 0]
            for comp_id in comp_ids:
                component_somas.setdefault(int(comp_id), set()).add(int(cid))

        d_ref = min(y_dim, x_dim) / 4.0
        edt = distance_transform_edt(~soma_bool).astype(np.float32)

        neighbors_26 = [(dz, dy, dx)
                        for dz in (-1, 0, 1)
                        for dy in (-1, 0, 1)
                        for dx in (-1, 0, 1)
                        if not (dz == 0 and dy == 0 and dx == 0)]

        def _make_slices(dz, dy, dx):
            src_s = [slice(None), slice(None), slice(None)]
            dst_s = [slice(None), slice(None), slice(None)]

            if dz > 0:
                src_s[0] = slice(0, z_dim - dz)
                dst_s[0] = slice(dz, z_dim)
            elif dz < 0:
                src_s[0] = slice(-dz, z_dim)
                dst_s[0] = slice(0, z_dim + dz)

            if dy > 0:
                src_s[1] = slice(0, y_dim - dy)
                dst_s[1] = slice(dy, y_dim)
            elif dy < 0:
                src_s[1] = slice(-dy, y_dim)
                dst_s[1] = slice(0, y_dim + dy)

            if dx > 0:
                src_s[2] = slice(0, x_dim - dx)
                dst_s[2] = slice(dx, x_dim)
            elif dx < 0:
                src_s[2] = slice(-dx, x_dim)
                dst_s[2] = slice(0, x_dim + dx)

            return tuple(src_s), tuple(dst_s)

        def _shifted_abs_diff(dz, dy, dx):
            src_s, dst_s = _make_slices(dz, dy, dx)
            shifted = np.zeros_like(astro)
            shifted[dst_s] = astro[src_s]
            return np.abs(shifted - astro).astype(np.float32)

        diff_max = np.zeros((z_dim, y_dim, x_dim), dtype=np.float32)
        for dz, dy, dx in neighbors_26:
            diff = _shifted_abs_diff(dz, dy, dx)
            diff_max = np.maximum(diff_max, diff)

        threshold = 0.5 * diff_max

        diff_big_count = np.zeros((z_dim, y_dim, x_dim), dtype=np.int32)
        for dz, dy, dx in neighbors_26:
            diff = _shifted_abs_diff(dz, dy, dx)
            big = diff > threshold
            diff_big_count += big.astype(np.int32)

        n_total = len(neighbors_26)
        n_big = diff_big_count
        n_small = n_total - diff_big_count
        n_small = np.maximum(n_small, 0)

        narrow_bonus = n_small.astype(np.float32) / float(n_total)

        base_weight = diff_max * astro
        base_weight = base_weight / (1.0 + edt / d_ref)

        weight = base_weight * (1.0 + narrow_bonus)

        lap_xyz = laplace(astro).astype(np.float32)
        dip = np.clip(lap_xyz, 0.0, None)
        eps_dip = 1e-6
        dip_norm = dip / (astro + eps_dip)

        weight = weight / (1.0 + dip_norm)
        dip_hard = dip_norm > dip_hard_threshold
        weight[dip_hard] = 0.0

        weight = np.clip(weight, 0.0, None).astype(np.float32)

        weight_in_graph = weight[graph_mask]
        if len(weight_in_graph) > 0:
            weight_median = float(np.median(weight_in_graph))
        else:
            weight_median = 0.0

        forbidden_mask = weight < (weight_median * 0.3)
        forbidden_mask = forbidden_mask & graph_mask
        forbidden_mask = forbidden_mask & (~soma_bool)

        eps = 1e-6
        z_step_penalty = 2.0

        final_cells = np.zeros_like(soma_mask, dtype=np.int32)

        for cid in unique_ids:
            if cid in soma_coords_by_cid:
                coords = soma_coords_by_cid[cid]
                for c in coords:
                    final_cells[int(c[0]), int(c[1]), int(c[2])] = int(cid)

        for comp_id in range(1, num_comp + 1):
            comp_somas = component_somas.get(comp_id, set())
            if len(comp_somas) == 0:
                continue

            if len(comp_somas) == 1:
                cid = next(iter(comp_somas))
                comp_mask = (labeled == comp_id)
                assign_mask = comp_mask & (~forbidden_mask)
                final_cells[assign_mask] = int(cid)
                continue

            comp_mask = (labeled == comp_id)
            comp_coords = np.argwhere(comp_mask)

            coord_to_idx = {}
            coords_list = []
            for i, c in enumerate(comp_coords):
                zz, yy, xx = int(c[0]), int(c[1]), int(c[2])
                coord_to_idx[(zz, yy, xx)] = i
                coords_list.append((zz, yy, xx))

            n_active = len(coords_list)
            dist_active = np.full(n_active, np.inf, dtype=np.float32)
            source_active = np.zeros(n_active, dtype=np.int32)

            heap = []
            for i, (zz, yy, xx) in enumerate(coords_list):
                if soma_bool[zz, yy, xx]:
                    dist_active[i] = 0.0
                    source_active[i] = int(soma_mask[zz, yy, xx])
                    heapq.heappush(heap, (0.0, i))

            while heap:
                d, i = heapq.heappop(heap)
                if d > dist_active[i]:
                    continue

                z, y, x = coords_list[i]
                w_here = weight[z, y, x]
                src_here = source_active[i]

                for dz, dy, dx in neighbors_26:
                    nz, ny, nx = z + dz, y + dy, x + dx
                    j = coord_to_idx.get((nz, ny, nx))
                    if j is None:
                        continue

                    if forbidden_mask[nz, ny, nx]:
                        continue

                    w_nb = weight[nz, ny, nx]
                    step_cost = 1.0 / (w_here * w_nb + eps)
                    if dz != 0:
                        step_cost *= z_step_penalty

                    nd = d + step_cost

                    if nd < dist_active[j]:
                        dist_active[j] = nd
                        source_active[j] = src_here
                        heapq.heappush(heap, (nd, j))

            for i, (zz, yy, xx) in enumerate(coords_list):
                src_id = int(source_active[i])
                if src_id > 0:
                    final_cells[zz, yy, xx] = src_id

        return final_cells

    def _segment_whole_cell_from_soma(self, astrocyte_stack, soma_mask):
        self.ctx.log(self.ctx.t("log.cells_title"))

        z_dim, y_dim, x_dim = astrocyte_stack.shape
        astro = astrocyte_stack.astype(np.float32)

        unique_ids = np.unique(soma_mask)[1:]
        if len(unique_ids) == 0:
            self.ctx.log(self.ctx.t("log.cells_no_soma"))
            self.ctx._last_astro_norm = None
            return np.zeros_like(astrocyte_stack, dtype=np.int32)

        image_mode = getattr(self.ctx, "image_mode", "confocal")

        if image_mode == "epifluorescence":
            astro_norm = self._preprocess_epi_log(astro, soma_mask)
            self.ctx._last_astro_norm = astro_norm
            return self._grow_cells_from_somas(astro_norm,
                                               soma_mask,
                                               unique_ids).astype(np.int32)

        median_all = float(np.median(astro[astro > 0])) if np.any(astro > 0) else 0.0
        if median_all > 0:
            eps = 1e-6
            k_max = 3.0
            for z in range(z_dim):
                slice_data = astro[z]
                nz = slice_data[slice_data > 0]
                if nz.size > 0:
                    med_z = float(np.median(nz))
                else:
                    med_z = 0.0
                med_z = max(med_z, eps)
                coef = median_all / med_z
                if coef > k_max:
                    coef = k_max
                astro[z] = slice_data * coef

        sigmas = [0.5, 1.0, 2.0, 3.0]

        structure_maps = []

        for sigma in sigmas:
            smoothed = gaussian_filter(astro, sigma=sigma)

            grad_z = sobel(smoothed, axis=0)
            grad_y = sobel(smoothed, axis=1)
            grad_x = sobel(smoothed, axis=2)

            grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)

            if np.max(grad_mag) > 0:
                grad_norm = grad_mag / np.max(grad_mag)
            else:
                grad_norm = grad_mag

            Gxx = gaussian_filter(grad_x**2, sigma=0.5)
            Gyy = gaussian_filter(grad_y**2, sigma=0.5)
            Gzz = gaussian_filter(grad_z**2, sigma=0.5)
            Gxy = gaussian_filter(grad_x * grad_y, sigma=0.5)
            Gxz = gaussian_filter(grad_x * grad_z, sigma=0.5)
            Gyz = gaussian_filter(grad_y * grad_z, sigma=0.5)

            trace = Gxx + Gyy + Gzz

            if np.max(trace) > 0:
                trace_norm = trace / np.max(trace)
            else:
                trace_norm = trace

            structure_map = grad_norm * trace_norm
            structure_maps.append(structure_map)

        max_structure = np.zeros_like(astro)
        min_structure = np.ones_like(astro) * 999
        mean_structure = np.zeros_like(astro)
        std_structure = np.zeros_like(astro)

        for i, smap in enumerate(structure_maps):
            max_structure = np.maximum(max_structure, smap)
            min_structure = np.minimum(min_structure, smap)
            mean_structure += smap

        mean_structure = mean_structure / len(structure_maps)

        for smap in structure_maps:
            std_structure += (smap - mean_structure)**2
        std_structure = np.sqrt(std_structure / len(structure_maps))

        grain_index = max_structure

        if np.max(std_structure) > 0:
            std_norm = std_structure / np.max(std_structure)
        else:
            std_norm = std_structure

        penalty = np.clip(std_norm * 0.3, 0, 0.3)
        grain_index = grain_index * (1 - penalty)

        small_scale = structure_maps[0]
        large_scale = structure_maps[-1]

        scale_ratio = large_scale / (small_scale + 1e-6)

        scale_penalty = np.clip((1 - scale_ratio) * 0.2, 0, 0.2)
        scale_penalty[scale_ratio > 0.5] = 0

        grain_index = grain_index * (1 - scale_penalty)

        grain_index = np.clip(grain_index, 0, 1)

        bg_grain = uniform_filter(grain_index, size=(3, 15, 15), mode='reflect')

        ratio = grain_index / (bg_grain + 1e-6)

        weight_map = 1.0 / (1.0 + np.exp(3.0 * (ratio - 1.0)))
        weight_map = np.clip(weight_map, 0.15, 1.0)

        weight_map = gaussian_filter(weight_map, sigma=(0.3, 0.5, 0.5))

        mip = np.max(astro, axis=0)
        sigma_bg = max(30, min(y_dim, x_dim) // 15)
        mip_bg = gaussian_filter(mip, sigma=sigma_bg)

        astro_clean = np.zeros_like(astro)
        for z in range(z_dim):
            base_clean = astro[z] - mip_bg * 0.7
            slice_weight = weight_map[z]
            astro_clean[z] = np.clip(base_clean * slice_weight, 0, None)

        uniform_mask = ratio <= 0.8

        for z in range(z_dim):
            slice_data = astro_clean[z]
            if np.max(slice_data) > 0:
                slice_bg = gaussian_filter(slice_data, sigma=max(20, min(y_dim, x_dim) // 20))
                mask_z = uniform_mask[z]
                if np.any(mask_z):
                    slice_corrected = slice_data.copy()
                    slice_corrected[mask_z] = np.clip(slice_data[mask_z] - slice_bg[mask_z] * 0.7, 0, None)
                    astro_clean[z] = slice_corrected

        astro_clean = median_filter(astro_clean, size=(1, 2, 2))

        if np.max(astro_clean) == 0:
            self.ctx.log(self.ctx.t("log.cells_no_signal"))
            self.ctx._last_astro_norm = None
            return soma_mask.astype(np.int32)

        p2, p98 = np.percentile(astro_clean[astro_clean > 0], [2, 98])
        astro_norm = np.clip((astro_clean - p2) / (p98 - p2 + 1e-6), 0, 1)

        self.ctx._last_astro_norm = astro_norm

        return self._grow_cells_from_somas(astro_norm,
                                           soma_mask,
                                           unique_ids).astype(np.int32)

    def _grow_cells_from_somas(self, astro_norm, soma_mask, unique_ids):
        z_dim, y_dim, x_dim = astro_norm.shape

        sigma_xy_bg = max(30, min(y_dim, x_dim) // 15)
        local_bg = uniform_filter(
            astro_norm.astype(np.float32),
            size=(3, sigma_xy_bg, sigma_xy_bg),
            mode='reflect'
        )

        local_bg = np.maximum(local_bg, 1e-6)

        contrast = astro_norm.astype(np.float32) / local_bg

        growth_mask = contrast > 1.2
        growth_mask = growth_mask & ~(soma_mask > 0)

        growth_mask = remove_small_objects(growth_mask, min_size=2)

        if growth_mask.max() == 0:
            self.ctx.log(self.ctx.t("log.cells_no_growth"))
            return soma_mask.astype(np.int32)

        final_cells = self._grow_cells_by_components(
            astro_norm, soma_mask, growth_mask, unique_ids
        )

        cleaned = np.zeros_like(final_cells)
        props = measure.regionprops(final_cells)

        min_cell_size = max(30, int(self.ctx.min_cell_volume))

        for prop in props:
            cid = prop.label
            if cid in unique_ids and prop.area > min_cell_size:
                cleaned[final_cells == cid] = cid

        final_count = len(np.unique(cleaned)) - 1
        self.ctx.log(self.ctx.t("log.cells_found", n=final_count))

        return cleaned.astype(np.int32)

    def _preprocess_epi_log(self, astro, soma_mask):
        from scipy.ndimage import laplace

        sigma_xy = 1.0
        sigma_z = 1.0

        blurred = gaussian_filter(
            astro,
            sigma=(sigma_z, sigma_xy, sigma_xy),
            mode='reflect'
        )

        log_map = -laplace(blurred)
        log_pos = np.clip(log_map, 0, None)

        nz = log_pos[log_pos > 0]
        if len(nz) == 0:
            return np.zeros_like(astro, dtype=np.float32)

        p2, p98 = np.percentile(nz, [2, 98])
        astro_norm = np.clip((log_pos - p2) / (p98 - p2 + 1e-6), 0, 1)
        astro_norm = self._apply_gamma(astro_norm, gamma=2.0)

        z_dim, y_dim, x_dim = astro_norm.shape
        sigma_s = 1.5
        k = 1.5
        sigma_l = max(15, min(y_dim, x_dim) // 15)

        lap_coh = gaussian_filter(
            log_map,
            sigma=(0.0, sigma_s, sigma_s),
            mode='reflect'
        )
        base_level = gaussian_filter(
            np.abs(lap_coh),
            sigma=(0.0, sigma_l, sigma_l),
            mode='reflect'
        )
        mask = np.abs(lap_coh) > (k * base_level)

        soma_bool = (soma_mask > 0)
        mask = mask | soma_bool

        astro_disc = astro_norm * mask.astype(np.float32)

        nz2 = astro_disc[astro_disc > 0]
        if len(nz2) == 0:
            return np.zeros_like(astro, dtype=np.float32)

        p2b, p98b = np.percentile(nz2, [2, 98])
        astro_disc = np.clip((astro_disc - p2b) / (p98b - p2b + 1e-6), 0, 1)
        return astro_disc.astype(np.float32)

    def _apply_gamma(self, astro_norm, gamma=2.0):
        if gamma == 1.0:
            return astro_norm.astype(np.float32)
        corrected = np.power(np.clip(astro_norm, 0.0, 1.0), gamma)
        return corrected.astype(np.float32)

class PostFiltering:
    def __init__(self, ctx, segmenter):
        self.ctx = ctx
        self.segmenter = segmenter

    def manual_cell_postfiltering(self, cell_mask, soma_mask, astro_norm, base_name, out_dir):
        self.ctx.log(self.ctx.t("log.filter_title"))

        z_dim, y_dim, x_dim = astro_norm.shape

        p0, p99 = np.percentile(astro_norm[astro_norm > 0], [1, 99])
        astro_norm_disp = np.clip((astro_norm - p0) / (p99 - p0 + 1e-6), 0, 1)

        mip = np.max(astro_norm_disp, axis=0)

        unique_cells = np.unique(cell_mask)[1:]
        if len(unique_cells) == 0:
            self.ctx.log(self.ctx.t("log.filter_no_cells"))
            return cell_mask

        cell_state = {int(cid): 0 for cid in unique_cells}

        original_cell_masks = {}
        original_soma_masks = {}
        for cid in unique_cells:
            cid_int = int(cid)
            original_cell_masks[cid_int] = (cell_mask == cid_int)
            original_soma_masks[cid_int] = (soma_mask == cid_int)

        cell_mip = np.zeros((y_dim, x_dim), dtype=np.int32)
        for z in range(z_dim):
            cell_mip[cell_mask[z] > 0] = cell_mask[z][cell_mask[z] > 0]

        cell_colors = {}
        unique_cells_list = list(unique_cells)
        for i, cid in enumerate(unique_cells_list):
            hue = (i * 0.618033988749895) % 1.0
            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
            cell_colors[int(cid)] = (r, g, b)

        res = self.ctx.bridge._ask(
            "cell_filter",
            mip=mip,
            cell_mip=cell_mip,
            colors=cell_colors,
        )
        if res is None or res.get("skip"):
            self.ctx.log(self.ctx.t("log.filter_skipped"))
            return cell_mask

        selected_cells = res.get("states", {})
        if not selected_cells:
            self.ctx.log(self.ctx.t("log.filter_no_selection"))
            return cell_mask

        for cid, state in selected_cells.items():
            cell_state[int(cid)] = state

        self.ctx.log(self.ctx.t("log.filter_resegment_n", n=len(selected_cells)))

        new_cell_mask = cell_mask.copy()
        resegmented_masks = {}

        for cid, problem_type in selected_cells.items():
            cell_old_mask = original_cell_masks[cid]
            soma_for_cell = original_soma_masks[cid]

            if not np.any(soma_for_cell):
                resegmented_masks[cid] = cell_old_mask
                continue

            coords = np.argwhere(cell_old_mask)
            if len(coords) == 0:
                resegmented_masks[cid] = cell_old_mask
                continue

            z_min = max(0, coords[:, 0].min() - 15)
            z_max = min(z_dim, coords[:, 0].max() + 15)
            y_min = max(0, coords[:, 1].min() - 60)
            y_max = min(y_dim, coords[:, 1].max() + 60)
            x_min = max(0, coords[:, 2].min() - 60)
            x_max = min(x_dim, coords[:, 2].max() + 60)

            local_astro_norm = astro_norm[z_min:z_max, y_min:y_max, x_min:x_max]
            local_soma = soma_for_cell[z_min:z_max, y_min:y_max, x_min:x_max]

            if problem_type == 1:
                mode = "strict"
            else:
                mode = "soft"

            new_local = self._resegment_local(
                local_astro_norm, local_soma, cid, mode
            )

            if np.any(new_local):
                new_cell_mask[cell_old_mask] = 0
                new_cell_mask[z_min:z_max, y_min:y_max, x_min:x_max][new_local > 0] = cid
                resegmented_masks[cid] = (new_cell_mask == cid)
            else:
                resegmented_masks[cid] = cell_old_mask

        self.ctx.log(self.ctx.t("log.filter_review"))

        new_masks = {cid: m for cid, m in resegmented_masks.items()}
        old_masks = {cid: original_cell_masks[cid] for cid in selected_cells}
        final_decision = self.ctx.bridge._ask(
            "resegment_review",
            mip=mip,
            new_masks=new_masks,
            old_masks=old_masks,
        )
        if final_decision is None:
            final_decision = {cid: None for cid in selected_cells}

        final_cell_mask = cell_mask.copy()
        kept_count = 0
        reverted_count = 0
        deleted_count = 0

        for cid, decision in final_decision.items():
            if decision == 'keep':
                kept_count += 1
            elif decision == 'delete':
                final_cell_mask[original_cell_masks[cid]] = 0
                deleted_count += 1
            else:
                final_cell_mask[original_cell_masks[cid]] = cid
                reverted_count += 1

        self.ctx.log(self.ctx.t("log.filter_kept", n=kept_count))
        self.ctx.log(self.ctx.t("log.filter_reverted", n=reverted_count))
        self.ctx.log(self.ctx.t("log.filter_deleted", n=deleted_count))

        decisions_dir = os.path.join(out_dir, "manual_filtering")
        os.makedirs(decisions_dir, exist_ok=True)

        decisions_file = os.path.join(decisions_dir, f"{base_name}_manual_filtering.json")
        filter_data = {
            'timestamp': time.time(),
            'cells': {
                str(cid): {
                    'problem': 'overgrown' if cell_state[cid] == 1 else 'undergrown',
                    'decision': final_decision[cid]
                }
                for cid in selected_cells.keys()
            },
            'stats': {
                'kept': kept_count,
                'reverted': reverted_count,
                'deleted': deleted_count
            }
        }
        with open(decisions_file, 'w') as f:
            json.dump(filter_data, f, indent=2)

        return final_cell_mask

    def _resegment_local(self, local_astro_norm, local_soma, cell_id, mode):
        from scipy.ndimage import binary_erosion

        z_dim, y_dim, x_dim = local_astro_norm.shape

        if z_dim <= 1 or y_dim <= 1 or x_dim <= 1:
            return np.zeros_like(local_soma, dtype=bool)

        astro_norm = local_astro_norm.astype(np.float32)

        nz = astro_norm[astro_norm > 0]
        if len(nz) == 0:
            return np.zeros_like(local_soma, dtype=bool)

        if mode == "strict":
            gamma = 3.0
            contrast_threshold = 2.0
            dip_hard_threshold = 0.3
            sigma_factor = 1.5
        else:
            gamma = 0.33
            contrast_threshold = 0.8
            dip_hard_threshold = 0.7
            sigma_factor = 0.7

        astro_norm_g = np.power(np.clip(astro_norm, 0.0, 1.0), gamma)

        sigma_xy_bg = max(15, int(min(y_dim, x_dim) // 15 * sigma_factor))
        local_bg = uniform_filter(
            astro_norm_g,
            size=(3, sigma_xy_bg, sigma_xy_bg),
            mode='reflect'
        )
        local_bg = np.maximum(local_bg, 1e-6)

        contrast = astro_norm_g / local_bg

        growth_mask = contrast > contrast_threshold
        growth_mask = growth_mask & ~(local_soma > 0)
        growth_mask = remove_small_objects(growth_mask, min_size=2)

        if growth_mask.max() == 0:
            return np.zeros_like(local_soma, dtype=bool)

        soma_int = local_soma.astype(np.int32)

        result = self.segmenter._grow_cells_by_components(
            astro_norm_g,
            soma_int,
            growth_mask,
            [int(cell_id)],
            dip_hard_threshold=dip_hard_threshold,
        )

        cell_mask = (result == int(cell_id))

        if not np.any(cell_mask):
            return np.zeros_like(local_soma, dtype=bool)

        if mode == "strict":
            eroded = binary_erosion(cell_mask, iterations=1)
            if np.any(eroded):
                cell_mask = eroded
        else:
            dilated = binary_dilation(cell_mask, iterations=1)
            dilated = dilated & growth_mask
            cell_mask = cell_mask | dilated

        return cell_mask

class Measurements:
    def __init__(self, ctx):
        self.ctx = ctx

    def _get_cell_suffix(self):
        return self.ctx.cell_type

    def _get_dim_suffix(self):
        return "2d" if self.ctx.is_2d else "3d"

    def _make_progress_logger(self, total, prefix=''):
        ctx = self.ctx

        class ProgressLogger:
            def __init__(self, total, prefix, verbose_each):
                self.total = total
                self.prefix = prefix
                self.n = 0
                self.start_time = time.time()
                self.last_log_time = 0.0
                self.last_logged_percent = -1
                self.verbose_each = verbose_each

            def _emit(self, msg):
                ctx.log(msg)

            def update(self, n=1):
                self.n += n
                now = time.time()
                elapsed = now - self.start_time

                percent = int(100 * self.n / self.total) if self.total > 0 else 100

                if self.n > 0 and self.total > 0:
                    eta = elapsed / self.n * (self.total - self.n)
                    eta_str = ctx.t("log.progress_eta", eta=eta)
                else:
                    eta_str = ""

                is_last = (self.n == self.total)
                time_ok = (now - self.last_log_time) >= 1.0
                percent_ok = (percent // 10) != (self.last_logged_percent // 10)

                if is_last or self.verbose_each or time_ok or percent_ok:
                    if is_last:
                        self._emit(ctx.t("log.progress_done",
                                         prefix=self.prefix,
                                         n=self.n,
                                         total=self.total,
                                         elapsed=elapsed))
                    else:
                        self._emit(ctx.t("log.progress_update",
                                         prefix=self.prefix,
                                         n=self.n,
                                         total=self.total,
                                         percent=percent,
                                         elapsed=elapsed,
                                         eta=eta_str))
                    self.last_log_time = now
                    self.last_logged_percent = percent

            def close(self):
                elapsed = time.time() - self.start_time
                if self.n < self.total:
                    self._emit(ctx.t("log.progress_closed_early",
                                     prefix=self.prefix,
                                     n=self.n,
                                     total=self.total,
                                     elapsed=elapsed))

        verbose_each = (total <= 20)
        return ProgressLogger(total, prefix, verbose_each)

    def measure_cells_2d(self, labels, intensity, nuclei, filename):
        if labels.shape != intensity.shape:
            self.ctx.log(self.ctx.t(
                "log.measurement_warn_shape",
                fn="measure_cells_2d",
                ls=labels.shape,
                is_=intensity.shape,
            ))

        if labels.max() == 0:
            return pd.DataFrame()

        self.ctx.log(self.ctx.t("log.measure_2d_title"))

        labels_2d = labels[0] if labels.ndim == 3 else labels
        intensity_2d = intensity[0] if intensity.ndim == 3 else intensity
        nuclei_2d = nuclei[0] if nuclei.ndim == 3 else nuclei

        props = measure.regionprops(labels_2d, intensity_image=intensity_2d)

        if not self.ctx.scale_initialized:
            self.ctx.xy_scale_um = 0.5
            self.ctx.z_scale_um = 2.0

        xy_um = self.ctx.xy_scale_um * self.ctx.downsample
        pixel_area_um2 = xy_um * xy_um

        data = []
        total = len(props)

        logger = self._make_progress_logger(total, prefix='    ')

        for idx, p in enumerate(props):
            cid = p.label
            cell_mask_2d = (labels_2d == cid)
            cy, cx = p.centroid

            cell_area_px = p.area
            cell_intensity = intensity_2d[cell_mask_2d]
            cell_mean_int = np.mean(cell_intensity) if len(cell_intensity) > 0 else 0
            cell_solidity = p.solidity if hasattr(p, 'solidity') else 1.0
            cell_perimeter = p.perimeter * xy_um if hasattr(p, 'perimeter') else 0

            if hasattr(p, 'inertia_tensor_eigenvalues'):
                ev = p.inertia_tensor_eigenvalues
                if len(ev) >= 2 and ev[0] > 0:
                    elongation = ev[-1] / ev[0]
                else:
                    elongation = 1
            else:
                elongation = 1

            if self.ctx._soma_labels.ndim == 3:
                soma_mask_2d = (self.ctx._soma_labels[0] == cid)
            else:
                soma_mask_2d = (self.ctx._soma_labels == cid)

            soma_area_px = np.sum(soma_mask_2d)
            if soma_area_px > 0:
                soma_intensity = intensity_2d[soma_mask_2d]
                soma_mean_int = np.mean(soma_intensity)
                soma_props = measure.regionprops(soma_mask_2d.astype(int), intensity_image=intensity_2d)
                soma_solidity = soma_props[0].solidity if len(soma_props) > 0 else 1.0
                soma_perimeter = soma_props[0].perimeter * xy_um if len(soma_props) > 0 else 0
            else:
                soma_mean_int = 0
                soma_solidity = 0
                soma_perimeter = 0

            nuc_area_px = 0
            nuc_mean_int = 0
            nuc_solidity = 0
            nuc_perimeter = 0

            if np.any(nuclei_2d > 0):
                nuclei_in_cell = np.unique(nuclei_2d[cell_mask_2d])
                nuclei_in_cell = nuclei_in_cell[nuclei_in_cell > 0]
                if len(nuclei_in_cell) > 0:
                    nuc_mask_2d = (nuclei_2d == nuclei_in_cell[0])
                    nuc_area_px = np.sum(nuc_mask_2d)
                    nuc_intensity = intensity_2d[nuc_mask_2d]
                    nuc_mean_int = np.mean(nuc_intensity)
                    nuc_props = measure.regionprops(nuc_mask_2d.astype(int), intensity_image=intensity_2d)
                    nuc_solidity = nuc_props[0].solidity if len(nuc_props) > 0 else 1.0
                    nuc_perimeter = nuc_props[0].perimeter * xy_um if len(nuc_props) > 0 else 0

            process_count = 0
            process_max_radius = 0
            if soma_area_px > 0:
                res = self.count_processes_ellipsoid_method(
                    labels, self.ctx._soma_labels
                )
                process_count = res['process_max_count']
                process_max_radius = res['process_max_radius']

            cell_area_um2 = cell_area_px * pixel_area_um2
            nuc_area_um2 = nuc_area_px * pixel_area_um2
            cyto_area_um2 = cell_area_um2 - nuc_area_um2
            if cyto_area_um2 > 0:
                nuc_cyto_ratio = nuc_area_um2 / cyto_area_um2
            else:
                nuc_cyto_ratio = 0

            data.append({
                'cell_id': cid,
                'filename': os.path.basename(filename),
                'centroid_y': cy,
                'centroid_x': cx,
                'cell_area_um2': cell_area_um2,
                'soma_area_um2': soma_area_px * pixel_area_um2,
                'nucleus_area_um2': nuc_area_um2,
                'cell_solidity': cell_solidity,
                'soma_solidity': soma_solidity,
                'nucleus_solidity': nuc_solidity,
                'elongation': elongation,
                'cell_mean_intensity': cell_mean_int,
                'soma_mean_intensity': soma_mean_int,
                'nucleus_mean_intensity': nuc_mean_int,
                'process_max_count': process_count,
                'process_max_radius': process_max_radius,
                'nuclear_cytoplasmic_ratio': nuc_cyto_ratio,
                'cell_perimeter_um': cell_perimeter,
                'soma_perimeter_um': soma_perimeter,
                'nucleus_perimeter_um': nuc_perimeter,
            })

            logger.update(1)

        logger.close()
        return pd.DataFrame(data)

    def measure_cells_3d(self, labels, intensity, nuclei, filename):
        if labels.shape != intensity.shape:
            self.ctx.log(self.ctx.t(
                "log.measurement_warn_shape",
                fn="measure_cells_3d",
                ls=labels.shape,
                is_=intensity.shape,
            ))

        if self.ctx.is_2d:
            return self.measure_cells_2d(labels, intensity, nuclei, filename)

        if labels.max() == 0:
            return pd.DataFrame()

        self.ctx.log(self.ctx.t("log.measure_3d_title"))

        props = measure.regionprops(labels, intensity_image=intensity)

        if not self.ctx.scale_initialized:
            self.ctx.xy_scale_um = 0.5
            self.ctx.z_scale_um = 2.0

        xy_um = self.ctx.xy_scale_um * self.ctx.downsample
        z_um = self.ctx.z_scale_um
        voxel_vol_um3 = xy_um * xy_um * z_um

        data = []
        total = len(props)

        logger = self._make_progress_logger(total, prefix='    ')

        for idx, p in enumerate(props):
            cid = p.label
            cell_mask = (labels == cid)
            cz, cy, cx = p.centroid

            cell_vol_vox = p.area
            cell_intensity = intensity[cell_mask]
            cell_mean_int = np.mean(cell_intensity) if len(cell_intensity) > 0 else 0
            cell_solidity = p.solidity if hasattr(p, 'solidity') else 1.0

            if hasattr(p, 'inertia_tensor_eigenvalues'):
                ev = p.inertia_tensor_eigenvalues
                if len(ev) >= 3 and ev[0] > 0:
                    elongation = ev[2] / ev[0]
                else:
                    elongation = 1
            else:
                elongation = 1

            cell_mip = np.max(cell_mask.astype(np.uint8), axis=0)
            cell_mip_perim = measure.perimeter(cell_mip) * xy_um

            soma_mask = (self.ctx._soma_labels == cid)
            soma_vol_vox = np.sum(soma_mask)
            if soma_vol_vox > 0:
                soma_intensity = intensity[soma_mask]
                soma_mean_int = np.mean(soma_intensity)
                soma_props = measure.regionprops(soma_mask.astype(int), intensity_image=intensity)
                soma_solidity = soma_props[0].solidity if len(soma_props) > 0 else 1.0
                soma_mip = np.max(soma_mask.astype(np.uint8), axis=0)
                soma_mip_perim = measure.perimeter(soma_mip) * xy_um
            else:
                soma_mean_int = 0
                soma_solidity = 0
                soma_mip_perim = 0

            nuc_vol_vox = 0
            nuc_mean_int = 0
            nuc_solidity = 0
            nuc_mip_perim = 0

            if np.any(nuclei > 0):
                nuclei_in_cell = np.unique(nuclei[cell_mask])
                nuclei_in_cell = nuclei_in_cell[nuclei_in_cell > 0]
                if len(nuclei_in_cell) > 0:
                    nuc_mask = (nuclei == nuclei_in_cell[0])
                    nuc_vol_vox = np.sum(nuc_mask)
                    nuc_intensity = intensity[nuc_mask]
                    nuc_mean_int = np.mean(nuc_intensity)
                    nuc_props = measure.regionprops(nuc_mask.astype(int), intensity_image=intensity)
                    nuc_solidity = nuc_props[0].solidity if len(nuc_props) > 0 else 1.0
                    nuc_mip = np.max(nuc_mask.astype(np.uint8), axis=0)
                    nuc_mip_perim = measure.perimeter(nuc_mip) * xy_um

            process_count = 0
            process_max_radius = 0
            if soma_vol_vox > 0:
                res = self.count_processes_ellipsoid_method(cell_mask, soma_mask)
                process_count = res['process_max_count']
                process_max_radius = res['process_max_radius']

            cell_volume_um3 = cell_vol_vox * voxel_vol_um3
            nuc_volume_um3 = nuc_vol_vox * voxel_vol_um3
            cyto_volume_um3 = cell_volume_um3 - nuc_volume_um3
            if cyto_volume_um3 > 0:
                nuc_cyto_ratio = nuc_volume_um3 / cyto_volume_um3
            else:
                nuc_cyto_ratio = 0

            data.append({
                'cell_id': cid,
                'filename': os.path.basename(filename),
                'centroid_z': cz,
                'centroid_y': cy,
                'centroid_x': cx,
                'cell_volume_um3': cell_volume_um3,
                'soma_volume_um3': soma_vol_vox * voxel_vol_um3,
                'nucleus_volume_um3': nuc_volume_um3,
                'cell_solidity': cell_solidity,
                'soma_solidity': soma_solidity,
                'nucleus_solidity': nuc_solidity,
                'elongation': elongation,
                'cell_mean_intensity': cell_mean_int,
                'soma_mean_intensity': soma_mean_int,
                'nucleus_mean_intensity': nuc_mean_int,
                'process_max_count': process_count,
                'process_max_radius': process_max_radius,
                'nuclear_cytoplasmic_ratio': nuc_cyto_ratio,
                'cell_MIP_projection_perimeter_um': cell_mip_perim,
                'soma_MIP_projection_perimeter_um': soma_mip_perim,
                'nucleus_MIP_projection_perimeter_um': nuc_mip_perim,
            })

            logger.update(1)

        logger.close()
        return pd.DataFrame(data)

    def count_processes_ellipsoid_method(self, cell_mask, nucleus_mask):
        nucleus_props = measure.regionprops(nucleus_mask.astype(int))[0]
        nz, ny, nx = nucleus_props.centroid
        nz, ny, nx = int(round(nz)), int(round(ny)), int(round(nx))

        nucleus_coords = np.argwhere(nucleus_mask)
        if len(nucleus_coords) > 0:
            z_range = nucleus_coords[:, 0].max() - nucleus_coords[:, 0].min() + 1
            y_range = nucleus_coords[:, 1].max() - nucleus_coords[:, 1].min() + 1
            x_range = nucleus_coords[:, 2].max() - nucleus_coords[:, 2].min() + 1

            soma_radius_z = max(3, z_range // 2 + 2)
            soma_radius_y = max(5, y_range // 2 + 3)
            soma_radius_x = max(5, x_range // 2 + 3)
        else:
            soma_radius_z, soma_radius_y, soma_radius_x = 3, 5, 5

        if self.ctx.scale_initialized:
            xy_scale = self.ctx.xy_scale_um * self.ctx.downsample
            z_scale = self.ctx.z_scale_um

            max_radius_um = 50
            max_radius_px = int(max_radius_um / xy_scale)
            min_radius_px = int(10 / xy_scale)

            step_um = 2.5
            step_px = max(2, int(step_um / xy_scale))

            radii_xy = np.arange(min_radius_px, max_radius_px + step_px, step_px)
        else:
            z_scale = 2.0
            radii_xy = np.arange(5, 40, 2)

        process_counts = []
        z_dim, y_dim, x_dim = cell_mask.shape

        for r_xy in radii_xy:
            if self.ctx.is_2d:
                r_z = 0
                ellipsoid = np.zeros((z_dim, y_dim, x_dim), dtype=bool)
                y_min = max(0, ny - r_xy)
                y_max = min(y_dim, ny + r_xy + 1)
                x_min = max(0, nx - r_xy)
                x_max = min(x_dim, nx + r_xy + 1)
                for y in range(y_min, y_max):
                    dy = (y - ny) / r_xy
                    dy2 = dy * dy
                    if dy2 > 1:
                        continue
                    for x in range(x_min, x_max):
                        dx = (x - nx) / r_xy
                        if dx * dx + dy2 <= 1:
                            ellipsoid[0, y, x] = True
            else:
                r_z = max(2, int(r_xy / z_scale))
                ellipsoid = np.zeros((z_dim, y_dim, x_dim), dtype=bool)
                z_min = max(0, nz - r_z)
                z_max = min(z_dim, nz + r_z + 1)
                y_min = max(0, ny - r_xy)
                y_max = min(y_dim, ny + r_xy + 1)
                x_min = max(0, nx - r_xy)
                x_max = min(x_dim, nx + r_xy + 1)
                for z in range(z_min, z_max):
                    dz = (z - nz) / r_z
                    dz2 = dz * dz
                    if dz2 > 1:
                        continue
                    for y in range(y_min, y_max):
                        dy = (y - ny) / r_xy
                        dy2 = dy * dy
                        if dz2 + dy2 > 1:
                            continue
                        for x in range(x_min, x_max):
                            dx = (x - nx) / r_xy
                            if dx * dx + dy2 + dz2 <= 1:
                                ellipsoid[z, y, x] = True

            soma_mask = binary_dilation(nucleus_mask, iterations=2)
            intersection = ellipsoid & cell_mask & ~soma_mask

            if intersection.sum() > 5:
                labeled_intersect, num = label(intersection, structure=np.ones((3, 3, 3)))
                props = measure.regionprops(labeled_intersect)
                valid_components = 0
                component_sizes = []
                for prop in props:
                    if prop.area >= 5:
                        if hasattr(prop, 'inertia_tensor_eigenvalues'):
                            ev = prop.inertia_tensor_eigenvalues
                            if len(ev) >= 3 and ev[0] > 0:
                                elongation = ev[2] / ev[0]
                                if elongation > 1.5:
                                    valid_components += 1
                                    component_sizes.append(prop.area)
                            else:
                                valid_components += 1
                                component_sizes.append(prop.area)
                        else:
                            valid_components += 1
                            component_sizes.append(prop.area)
                process_counts.append({
                    'count': valid_components,
                    'sizes': component_sizes,
                    'total_voxels': intersection.sum()
                })
            else:
                process_counts.append({'count': 0, 'sizes': [], 'total_voxels': 0})

        if len(process_counts) > 0:
            counts = [p['count'] for p in process_counts]
            sizes_flat = [s for p in process_counts for s in p['sizes']]
            max_processes = max(counts)
            max_idx = np.argmax(counts)
            max_radius = radii_xy[max_idx] if max_idx < len(radii_xy) else 0
            mid_idx = len(radii_xy) // 2
            if len(counts) > 5:
                avg_mid = np.mean(counts[mid_idx-2:mid_idx+3])
            else:
                avg_mid = np.mean(counts) if counts else 0
            integral = np.trapezoid(counts, dx=2) if len(counts) > 1 else 0
            if max_processes > 0 and len(counts) > 5:
                decay_ratio = counts[-1] / max_processes if max_processes > 0 else 0
            else:
                decay_ratio = 0
            if len(counts) > 3:
                derivatives = np.diff(counts)
                branching_complexity = np.std(derivatives) if len(derivatives) > 0 else 0
            else:
                branching_complexity = 0
            if max_processes > 0:
                half_max = max_processes * 0.5
                half_max_radius = 0
                for i, cnt in enumerate(counts):
                    if cnt < half_max:
                        half_max_radius = radii_xy[i] if i < len(radii_xy) else radii_xy[-1]
                        break
                if half_max_radius == 0:
                    half_max_radius = radii_xy[-1]
            else:
                half_max_radius = 0
            avg_process_size = np.mean(sizes_flat) if sizes_flat else 0
            total_process_voxels = sum([p['total_voxels'] for p in process_counts])
            cell_volume = np.sum(cell_mask)
            processes_to_cell_ratio = total_process_voxels / (cell_volume + 1)

            return {
                'process_max_count': max_processes,
                'process_max_radius': max_radius,
                'process_avg_mid': avg_mid,
                'process_integral': integral,
                'process_decay_ratio': decay_ratio,
                'process_branching_complexity': branching_complexity,
                'process_half_max_radius': half_max_radius,
                'process_avg_size': avg_process_size,
                'process_total_voxels': total_process_voxels,
                'process_to_cell_ratio': processes_to_cell_ratio,
                'process_profile': counts,
                'process_radii': radii_xy.tolist()
            }
        else:
            return {
                'process_max_count': 0, 'process_max_radius': 0, 'process_avg_mid': 0,
                'process_integral': 0, 'process_decay_ratio': 0, 'process_branching_complexity': 0,
                'process_half_max_radius': 0, 'process_avg_size': 0, 'process_total_voxels': 0,
                'process_to_cell_ratio': 0, 'process_profile': [], 'process_radii': []
            }

    def _append_to_all_measurements(self, measurements, filepath):
        if measurements is None or len(measurements) == 0:
            return

        suffix = self._get_cell_suffix()
        dim_suffix = self._get_dim_suffix()
        out_dir = os.path.join(os.path.dirname(filepath), f"all_results_{suffix}_{dim_suffix}")
        os.makedirs(out_dir, exist_ok=True)

        all_csv = os.path.join(out_dir, "all_measurements.csv")
        current_filename = os.path.basename(filepath)

        if os.path.exists(all_csv):
            try:
                existing = pd.read_csv(all_csv)
            except Exception as e:
                self.ctx.log(self.ctx.t("log.results_read_failed",
                                        path=all_csv, err=e))
                existing = pd.DataFrame()

            if len(existing) > 0 and 'filename' in existing.columns:
                existing = existing[existing['filename'] != current_filename]

            combined = pd.concat([existing, measurements], ignore_index=True)
        else:
            combined = measurements.copy()

        combined.to_csv(all_csv, index=False)
        self.ctx.log(self.ctx.t("log.results_appended",
                                path=all_csv, n=len(combined)))

class Visualization:
    def __init__(self, ctx):
        self.ctx = ctx

    def _save_cells_overlay_stack(self, astro, cells, base_name, out_dir):
        stacks_dir = os.path.join(out_dir, "3d_stacks")
        os.makedirs(stacks_dir, exist_ok=True)

        z_dim, y_dim, x_dim = astro.shape

        astro_norm = astro.astype(np.float32)
        nz = astro_norm[astro_norm > 0]
        if len(nz) > 0:
            p5, p95 = np.percentile(nz, [5, 95])
            if p95 - p5 < 1e-6:
                astro_norm = np.zeros_like(astro_norm)
            else:
                astro_norm = np.clip((astro_norm - p5) / (p95 - p5), 0, 1)
        else:
            astro_norm = np.zeros_like(astro_norm)

        rgb_stack = np.zeros((z_dim, y_dim, x_dim, 3), dtype=np.uint8)

        base_u8 = (astro_norm * 255).astype(np.uint8)
        for ch in range(3):
            rgb_stack[..., ch] = base_u8

        if cells.max() > 0:
            unique_cells = np.unique(cells)[1:]
            alpha = 0.65

            for i, cid in enumerate(unique_cells):
                hue = (i * 0.618033988749895) % 1.0
                r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 1.0)
                mask = (cells == cid)
                if not np.any(mask):
                    continue

                r_u8 = int(r * 255)
                g_u8 = int(g * 255)
                b_u8 = int(b * 255)

                for ch, col in enumerate((r_u8, g_u8, b_u8)):
                    channel = rgb_stack[..., ch].astype(np.float32)
                    channel[mask] = (
                            channel[mask] * (1.0 - alpha) + col * alpha
                    )
                    rgb_stack[..., ch] = channel.astype(np.uint8)

        tifffile.imwrite(
            os.path.join(stacks_dir, f"{base_name}_cells_overlay_rgb.tiff"),
            rgb_stack,
            compression='zlib'
        )

        self.ctx.log(self.ctx.t("log.viz_overlay_saved",
                                name=f"{base_name}_cells_overlay_rgb.tiff"))

    def save_3d_stacks(self, astro, dapi, cells, nuclei, base_name, out_dir):
        stacks_dir = os.path.join(out_dir, "3d_stacks")
        os.makedirs(stacks_dir, exist_ok=True)

        astro_norm = ((astro - astro.min()) / (astro.max() - astro.min() + 1e-6) * 65535).astype(np.uint16)
        tifffile.imwrite(
            os.path.join(stacks_dir, f"{base_name}_astro_3d.tiff"),
            astro_norm,
            compression='zlib'
        )

        dapi_norm = ((dapi - dapi.min()) / (dapi.max() - dapi.min() + 1e-6) * 65535).astype(np.uint16)
        tifffile.imwrite(
            os.path.join(stacks_dir, f"{base_name}_dapi_3d.tiff"),
            dapi_norm,
            compression='zlib'
        )

        if nuclei.max() > 0:
            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_nuclei_ids.tiff"),
                nuclei.astype(np.uint16),
                compression='zlib'
            )

        if cells.max() > 0:
            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_astrocytes_ids.tiff"),
                cells.astype(np.uint16),
                compression='zlib'
            )

            z_dim, y_dim, x_dim = cells.shape
            rgb_stack = np.zeros((z_dim, y_dim, x_dim, 3), dtype=np.uint8)

            unique_cells = np.unique(cells)[1:]
            for i, cell_id in enumerate(unique_cells):
                hue = (i * 0.618033988749895) % 1.0
                r, g, b = colorsys.hsv_to_rgb(hue, 0.8, 1.0)
                mask = cells == cell_id
                rgb_stack[mask] = [int(r*255), int(g*255), int(b*255)]

            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_astrocytes_rgb.tiff"),
                rgb_stack,
                compression='zlib'
            )

            self.ctx.log(self.ctx.t("log.viz_saved_n_cells",
                                    n=len(unique_cells)))

    def _save_slices_visualization(self, image_stack, mask, title, filename_suffix):
        out_dir = self.ctx._current_out_dir or '.'
        viz_dir = os.path.join(out_dir, "visualizations")
        os.makedirs(viz_dir, exist_ok=True)

        z_dim, y_dim, x_dim = image_stack.shape
        img_norm = image_stack.astype(np.float32)
        nz = img_norm[img_norm > 0]
        if len(nz) == 0:
            return
        p0, p99 = np.percentile(nz, [1, 99])
        img_norm = np.clip((img_norm - p0) / (p99 - p0 + 1e-6), 0, 1)

        unique_ids = np.unique(mask)[1:]
        if len(unique_ids) == 0:
            return

        z_occupied = np.where(np.any(mask, axis=(1, 2)))[0]
        if len(z_occupied) == 0:
            return

        nrows, ncols = 2, 5
        if len(z_occupied) > nrows * ncols:
            selected_indices = np.linspace(0, len(z_occupied) - 1, nrows * ncols, dtype=int)
            selected_z = z_occupied[selected_indices]
        else:
            selected_z = z_occupied

        cell_colors = {}
        for cid in unique_ids:
            hue = (cid * 0.618033988749895) % 1.0
            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
            cell_colors[cid] = (r, g, b)

        fig, axes = plt.subplots(nrows, ncols, figsize=(20, 8))
        axes = axes.flatten()

        for idx, z in enumerate(selected_z):
            ax = axes[idx]
            ax.imshow(img_norm[z], cmap='gray')

            for cid in unique_ids:
                mask_2d = (mask[z] == cid)
                if np.any(mask_2d):
                    r, g, b = cell_colors[cid]
                    overlay = np.zeros((y_dim, x_dim, 4))
                    overlay[mask_2d] = [r, g, b, 0.8]
                    ax.imshow(overlay)

            ax.set_title(self.ctx.t("log.slice_title", z=z), color='black')
            ax.axis('off')

        for i in range(len(selected_z), len(axes)):
            axes[i].axis('off')

        plt.suptitle(title, fontsize=16)
        plt.tight_layout()

        base_name = self.ctx._current_base_name or 'unknown'
        plt.savefig(
            os.path.join(viz_dir, f"{base_name}_{filename_suffix}.png"),
            dpi=200, bbox_inches='tight'
        )
        plt.close(fig)

    def save_visualizations(self, astro, dapi, cells, nuclei, base_name, out_dir):
        viz_dir = os.path.join(out_dir, "visualizations")
        os.makedirs(viz_dir, exist_ok=True)
        z_dim, y_dim, x_dim = astro.shape

        self.ctx.log(self.ctx.t("log.viz_final_start"))

        if cells.max() > 0:
            tifffile.imwrite(
                os.path.join(viz_dir, f"{base_name}_mask_astrocytes.tiff"),
                cells.astype(np.uint16),
                compression='zlib'
            )
            if self.ctx._soma_labels is not None:
                tifffile.imwrite(
                    os.path.join(viz_dir, f"{base_name}_mask_somas.tiff"),
                    self.ctx._soma_labels.astype(np.uint16),
                    compression='zlib'
                )
            if nuclei.max() > 0:
                tifffile.imwrite(
                    os.path.join(viz_dir, f"{base_name}_mask_nuclei.tiff"),
                    nuclei.astype(np.uint16),
                    compression='zlib'
                )
            rgb_stack = np.zeros((z_dim, y_dim, x_dim, 3), dtype=np.uint8)
            unique_cells = np.unique(cells)[1:]
            for i, cid in enumerate(unique_cells):
                hue = (i * 0.618033988749895) % 1.0
                r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 1.0)
                mask = cells == cid
                rgb_stack[mask] = [int(r*255), int(g*255), int(b*255)]
            tifffile.imwrite(
                os.path.join(viz_dir, f"{base_name}_mask_astrocytes_rgb.tiff"),
                rgb_stack,
                compression='zlib'
            )

        astro_norm = astro.astype(np.float32)
        nz = astro_norm[astro_norm > 0]
        if len(nz) > 0:
            p5, p95 = np.percentile(nz, [5, 95])
            astro_norm = np.clip((astro_norm - p5) / (p95 - p5 + 1e-6), 0, 1)
        astro_norm = np.power(astro_norm, 1.2)
        mip_bg = np.max(astro_norm, axis=0)

        if cells.max() > 0:
            bg_rgb = np.zeros((y_dim, x_dim, 3), dtype=np.uint8)
            bg_rgb[..., 0] = (mip_bg * 0.30).astype(np.uint8)
            bg_rgb[..., 1] = (mip_bg * 0.40).astype(np.uint8)
            bg_rgb[..., 2] = (mip_bg * 0.50).astype(np.uint8)

            overlay_rgba = np.zeros((y_dim, x_dim, 4), dtype=np.float32)
            overlay_rgba[..., :3] = bg_rgb / 255.0
            overlay_rgba[..., 3] = 1.0

            unique_cells = np.unique(cells)[1:]
            for i, cid in enumerate(unique_cells):
                hue = (i * 0.618033988749895) % 1.0
                r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 1.0)
                cell_mip_2d = np.max((cells == cid).astype(np.uint8), axis=0)

                for ch in range(3):
                    overlay_rgba[cell_mip_2d > 0, ch] = (overlay_rgba[cell_mip_2d > 0, ch] * 0.35) + (r * 0.65)
                overlay_rgba[cell_mip_2d > 0, 3] = 1.0

            final_mip_bg = (overlay_rgba[..., :3] * 255).astype(np.uint8)

            mip_rgb_black = np.zeros((y_dim, x_dim, 3), dtype=np.uint8)
            for i, cid in enumerate(unique_cells):
                hue = (i * 0.618033988749895) % 1.0
                r, g, b = colorsys.hsv_to_rgb(hue, 0.9, 1.0)
                cell_mip_2d = np.max((cells == cid).astype(np.uint8), axis=0)
                for ch in range(3):
                    mip_rgb_black[cell_mip_2d > 0, ch] = int([r, g, b][ch] * 220)

            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            ax.imshow(final_mip_bg)
            ax.axis('off')
            plt.savefig(os.path.join(viz_dir, f"{base_name}_mip_colored_bg.png"), bbox_inches='tight', dpi=300)
            plt.close()

            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            ax.imshow(mip_rgb_black)
            ax.axis('off')
            plt.savefig(os.path.join(viz_dir, f"{base_name}_mip_colored_black.png"), bbox_inches='tight', dpi=300)
            plt.close()

        self.ctx.log(self.ctx.t("log.viz_saved_to", path=viz_dir))

    def _save_deconvolved_stacks(self, astro, dapi, base_name, out_dir):
        deconv_dir = os.path.join(out_dir, "deconvolved")
        os.makedirs(deconv_dir, exist_ok=True)

        if astro is not None and astro.size > 0:
            astro_norm = astro.astype(np.float32)
            amax = np.max(astro_norm)
            if amax > 0:
                astro_u16 = (astro_norm / amax * 65535).astype(np.uint16)
            else:
                astro_u16 = astro_norm.astype(np.uint16)

            tifffile.imwrite(
                os.path.join(deconv_dir, f"{base_name}_cells_deconvolved.tiff"),
                astro_u16,
                compression='zlib'
            )
            self.ctx.log(self.ctx.t("log.viz_deconv_cells_saved",
                                    shape=astro.shape))

            try:
                mip = np.max(astro_norm, axis=0)
                nzm = mip[mip > 0]
                if len(nzm) > 0:
                    p0, p99 = np.percentile(nzm, [1, 99])
                    mip_norm = np.clip((mip - p0) / (p99 - p0 + 1e-6), 0, 1)
                else:
                    mip_norm = mip

                fig, ax = plt.subplots(1, 1, figsize=(8, 8))
                ax.imshow(mip_norm, cmap='gray')
                ax.set_title(self.ctx.t("log.deconv_cells_title", name=base_name), fontsize=12)
                ax.axis('off')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(deconv_dir, f"{base_name}_cells_deconvolved_mip.png"),
                    dpi=200, bbox_inches='tight'
                )
                plt.close(fig)
            except Exception as e:
                self.ctx.log(self.ctx.t("log.viz_mip_cells_failed", err=e))

        if dapi is not None and dapi.size > 0:
            dapi_norm = dapi.astype(np.float32)
            dmax = np.max(dapi_norm)
            if dmax > 0:
                dapi_u16 = (dapi_norm / dmax * 65535).astype(np.uint16)
            else:
                dapi_u16 = dapi_norm.astype(np.uint16)

            tifffile.imwrite(
                os.path.join(deconv_dir, f"{base_name}_dapi_deconvolved.tiff"),
                dapi_u16,
                compression='zlib'
            )
            self.ctx.log(self.ctx.t("log.viz_deconv_dapi_saved",
                                    shape=dapi.shape))

            try:
                mip = np.max(dapi_norm, axis=0)
                nzm = mip[mip > 0]
                if len(nzm) > 0:
                    p0, p99 = np.percentile(nzm, [1, 99])
                    mip_norm = np.clip((mip - p0) / (p99 - p0 + 1e-6), 0, 1)
                else:
                    mip_norm = mip

                fig, ax = plt.subplots(1, 1, figsize=(8, 8))
                ax.imshow(mip_norm, cmap='gray')
                ax.set_title(self.ctx.t("log.deconv_dapi_title", name=base_name), fontsize=12)
                ax.axis('off')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(deconv_dir, f"{base_name}_dapi_deconvolved_mip.png"),
                    dpi=200, bbox_inches='tight'
                )
                plt.close(fig)
            except Exception as e:
                self.ctx.log(self.ctx.t("log.viz_mip_dapi_failed", err=e))

class Analyzer:
    def __init__(self):
        self.ctx = AnalysisContext()
        self.file_io = FileIO(self.ctx)
        self.optics = Optics(self.ctx, self.file_io)
        self.deconv = Deconvolution(self.ctx)
        self.segmenter = Segmentation(self.ctx)
        self.postfilter = PostFiltering(self.ctx, self.segmenter)
        self.measurer = Measurements(self.ctx)
        self.viz = Visualization(self.ctx)

    def batch_process_files(self, files):
        self.ctx.log(self.ctx.t("log.batch_type_cells",
                                t=self.ctx.cell_type_ru))
        mode_ru = (self.ctx.t("mode.epifluorescence_full")
                   if self.ctx.image_mode == "epifluorescence"
                   else self.ctx.t("mode.confocal_full"))
        self.ctx.log(self.ctx.t("log.batch_image_mode", m=mode_ru))
        state_deconv = (self.ctx.t("state.on")
                        if self.ctx.deconvolve_images
                        else self.ctx.t("state.off"))
        self.ctx.log(self.ctx.t("log.batch_deconv", s=state_deconv))
        state_filter = (self.ctx.t("state.on")
                        if self.ctx._gui_do_filter
                        else self.ctx.t("state.off"))
        self.ctx.log(self.ctx.t("log.batch_filter", s=state_filter))

        if not files:
            return

        self.ctx.log(self.ctx.t("log.batch_files_count", n=len(files)))

        all_results = []
        start = time.time()

        for i, f in enumerate(files, 1):
            self.ctx.log(self.ctx.t("log.batch_file_n",
                                    i=i, total=len(files),
                                    name=os.path.basename(f)))
            m, _ = self.process_file(f)
            if m is not None and len(m) > 0:
                all_results.append(m)

        if all_results:
            combined = pd.concat(all_results, ignore_index=True)

            self.ctx.log(self.ctx.t("log.batch_files_done", n=len(files)))
            self.ctx.log(self.ctx.t("log.batch_cells_found", n=len(combined)))
            self.ctx.log(self.ctx.t("log.batch_total_time",
                                    t=time.time() - start))

            if 'nuclear_cytoplasmic_ratio' in combined.columns:
                self.ctx.log(self.ctx.t("log.batch_averages"))
                if self.ctx.is_2d:
                    if 'cell_area_um2' in combined.columns:
                        self.ctx.log(self.ctx.t("log.batch_avg_area",
                                                v=combined['cell_area_um2'].mean()))
                    if 'nucleus_area_um2' in combined.columns:
                        self.ctx.log(self.ctx.t("log.batch_avg_nuc_area",
                                                v=combined['nucleus_area_um2'].mean()))
                else:
                    if 'cell_volume_um3' in combined.columns:
                        self.ctx.log(self.ctx.t("log.batch_avg_volume",
                                                v=combined['cell_volume_um3'].mean()))
                    if 'nucleus_volume_um3' in combined.columns:
                        self.ctx.log(self.ctx.t("log.batch_avg_nuc_volume",
                                                v=combined['nucleus_volume_um3'].mean()))
                self.ctx.log(self.ctx.t("log.batch_avg_nc_ratio",
                                        v=combined['nuclear_cytoplasmic_ratio'].mean()))
                self.ctx.log(self.ctx.t("log.batch_avg_elong",
                                        v=combined['elongation'].mean()))
        else:
            self.ctx.log(self.ctx.t("log.batch_no_cells"))

    def batch_process(self):
        files = self.file_io.get_file_paths()
        if not files:
            return
        self.batch_process_files(files)

    def process_file(self, filepath):
        try:
            start_time = time.time()

            if not self.ctx.scale_initialized:
                self.file_io.extract_and_confirm_metadata(filepath)

            if not self.ctx._channels_initialized:
                self.optics.determine_channels(filepath)
                self.ctx._channels_initialized = True

            if self.ctx.deconvolve_images and not self.ctx._optics_initialized:
                self.optics.extract_and_confirm_optics(filepath)
                self.ctx._optics_initialized = True

            img = self.file_io.load_and_downsample(filepath)
            if img is None:
                return None, None

            astro = img[0]
            if self.ctx.has_dapi and img.shape[0] > 1:
                dapi = img[1]
            else:
                dapi = np.zeros_like(astro)
                self.ctx.log(self.ctx.t("log.dapi_not_used_msg"))

            if self.ctx.deconvolve_images:
                pixel_xy_um = self.ctx.xy_scale_um * self.ctx.downsample
                pixel_z_um = self.ctx.z_scale_um

                astro = self.deconv._deconvolve_3d_stack(astro, pixel_xy_um, pixel_z_um)

                if self.ctx.has_dapi and np.any(dapi > 0):
                    dapi = self.deconv._deconvolve_3d_stack_dapi(dapi, pixel_xy_um, pixel_z_um)

                self.ctx._deconvolution_applied = True
            else:
                self.ctx._deconvolution_applied = False

            base = os.path.splitext(os.path.basename(filepath))[0]
            suffix = self.measurer._get_cell_suffix()
            out_dir = os.path.join(os.path.dirname(filepath), f"{base}_{suffix}_analysis")
            os.makedirs(out_dir, exist_ok=True)

            self.ctx._current_base_name = base
            self.ctx._current_out_dir = out_dir

            self.viz._save_deconvolved_stacks(astro, dapi, base, out_dir)

            t2 = time.time()

            self.ctx.soma_segmented = self.segmenter.manual_soma_setting(astro)
            if self.ctx.soma_segmented.max() == 0:
                self.ctx.log(self.ctx.t("log.soma_not_marked"))
                return None, None
            self.viz._save_slices_visualization(astro, self.ctx.soma_segmented, self.ctx.t("log.title_soma_slices"), "soma_slices")
            self.viz.save_3d_stacks(astro, dapi, self.ctx.soma_segmented, np.zeros_like(dapi), base, out_dir)

            self.ctx.nuclei_segmented = self.segmenter.segment_nuclei_local_from_soma(self.ctx.soma_segmented, dapi)
            self.viz._save_slices_visualization(dapi, self.ctx.nuclei_segmented, self.ctx.t("log.title_nuclei_slices"), "nuclei_slices")
            self.viz.save_3d_stacks(astro, dapi, self.ctx.nuclei_segmented, np.zeros_like(dapi), base, out_dir)

            self.ctx.cells_segmented = self.segmenter._segment_whole_cell_from_soma(astro, self.ctx.soma_segmented)
            self.viz._save_slices_visualization(astro, self.ctx.cells_segmented, self.ctx.t("log.title_cells_slices"), "cells_slices")
            self.viz.save_3d_stacks(astro, dapi, self.ctx.cells_segmented, np.zeros_like(dapi), base, out_dir)

            self.ctx.log(self.ctx.t("log.segmentation_time",
                                    t=time.time() - t2))

            if self.ctx._gui_do_filter:
                astro_norm_for_filter = getattr(self.ctx, "_last_astro_norm", None)
                if astro_norm_for_filter is None:
                    self.ctx.log(self.ctx.t("log.filter_skip_no_norm"))
                else:
                    self.ctx.cells_segmented = self.postfilter.manual_cell_postfiltering(
                        self.ctx.cells_segmented,
                        self.ctx.soma_segmented,
                        astro_norm_for_filter,
                        base,
                        out_dir
                    )
            else:
                self.ctx.log(self.ctx.t("log.filter_skip"))

            self.viz._save_cells_overlay_stack(astro, self.ctx.cells_segmented, base, out_dir)

            t3 = time.time()
            measurements = self.measurer.measure_cells_3d(self.ctx.cells_segmented, astro, self.ctx.nuclei_segmented, filepath)
            self.ctx.log(self.ctx.t("log.measurement_time",
                                    t=time.time() - t3))

            self.viz.save_visualizations(astro, dapi, self.ctx.cells_segmented, self.ctx.nuclei_segmented, base, out_dir)

            if len(measurements) > 0:
                measurements.to_csv(
                    os.path.join(out_dir, f"{base}_measurements.csv"),
                    index=False
                )
                self.measurer._append_to_all_measurements(measurements, filepath)

            self.ctx.log(self.ctx.t("log.file_summary"))
            self.ctx.log(self.ctx.t("log.file_summary_cells", n=len(measurements)))
            self.ctx.log(self.ctx.t("log.file_summary_time",
                                    t=time.time() - start_time))

            return measurements, None

        except Exception as e:
            self.ctx.log(self.ctx.t("log.error_prefix", err=e))
            self.ctx.log(traceback.format_exc())
            return None, None

class GuiBridge:
    def __init__(self, app):
        self.app = app
        self._answer_q = queue.Queue()

    def log(self, text):
        self.app.append_log(str(text))

    def _ask(self, kind, **kwargs):
        self.app.root.after(0, lambda: self._dispatch(kind, kwargs))
        return self._answer_q.get()

    def _ok(self, value):
        self._answer_q.put(value)

    def _dispatch(self, kind, kwargs):
        try:
            if kind == "file_paths":         self._dialog_file_paths(kwargs)
            elif kind == "channel_selector": self._dialog_channel_selector(kwargs)
            elif kind == "dapi_channel":     self._dialog_dapi_channel(kwargs)
            elif kind == "soma_picker":      self._dialog_soma_picker(kwargs)
            elif kind == "soma_confirm":     self._dialog_soma_confirm(kwargs)
            elif kind == "slice_selector":   self._dialog_slice_selector(kwargs)
            elif kind == "cell_filter":      self._dialog_cell_filter(kwargs)
            elif kind == "resegment_review": self._dialog_resegment_review(kwargs)
            elif kind == "soma_repick":      self._dialog_soma_repick(kwargs)
            else:
                messagebox.showerror(self.app._t("app.title"),
                                     self.app._t("dlg.unknown", kind=kind))
                self._ok(None)
        except Exception as e:
            self.app.append_log(self.app._t("log.dialog_error", err=str(e)))
            self.app.append_log(traceback.format_exc())
            self._ok(None)

    @staticmethod
    def _to_uint8(img):
        img = np.asarray(img, dtype=np.float32)
        nz = img[img > 0]
        if nz.size == 0:
            return np.zeros(img.shape, dtype=np.uint8)
        p0, p99 = np.percentile(nz, [1, 99])
        if p99 - p0 < 1e-6:
            return np.zeros(img.shape, dtype=np.uint8)
        return np.clip((img - p0) / (p99 - p0) * 255.0, 0, 255).astype(np.uint8)

    @staticmethod
    def _pil_from_gray(arr_u8):
        return Image.fromarray(arr_u8, mode="L")

    def _fit_size(self, w, h, max_w=520, max_h=520):
        scale = min(1.0, max_w / max(w, 1), max_h / max(h, 1))
        return int(w * scale), int(h * scale), scale

    def _dialog_file_paths(self, kw):
        exts = kw.get("exts", ['.tif', '.tiff', '.nd2', '.oib', '.oif', '.oir'])
        panel = self.app.start_work_panel(self.app._t("dlg.files.title"))

        tk.Label(panel, text=self.app._t("dlg.files.one_or_dir"),
                 font=("", 11, "bold")).pack(anchor="w", pady=(0, 8))
        mode = tk.StringVar(value="one")
        tk.Radiobutton(panel, text=self.app._t("dlg.files.one"),
                       variable=mode, value="one").pack(anchor="w")
        tk.Radiobutton(panel, text=self.app._t("dlg.files.dir"),
                       variable=mode, value="dir").pack(anchor="w")

        selected = {"files": []}
        info_lbl = tk.Label(panel, text="", justify="left", wraplength=480, fg="#333")
        info_lbl.pack(anchor="w", pady=(8, 0))

        def choose():
            if mode.get() == "one":
                p = filedialog.askopenfilename(
                    filetypes=[
                        (self.app._t("dlg.files.images"),
                         "*.tif *.tiff *.nd2 *.czi *.lif *.ome.tif *.ome.tiff *.oib *.oif *.oir"),
                        ("TIFF", "*.tif *.tiff *.ome.tif *.ome.tiff"),
                        ("ND2",  "*.nd2"),
                        ("CZI",  "*.czi"),
                        ("LIF",  "*.lif"),
                        ("Olympus", "*.oib *.oif *.oir"),
                        (self.app._t("dlg.files.all"), "*.*"),
                    ])
                if not p:
                    return
                selected["files"] = [p]
            else:
                d = filedialog.askdirectory()
                if not d:
                    return
                selected["files"] = sorted(
                    os.path.join(d, f) for f in os.listdir(d)
                    if any(f.lower().endswith(e) for e in exts))
            info_lbl.config(text=self.app._t("dlg.files.selected_n",
                                             n=len(selected["files"])))

        tk.Button(panel, text=self.app._t("btn.choose"),
                  command=choose).pack(anchor="w", pady=(12, 0))

        def finish():
            if not selected["files"]:
                messagebox.showwarning(self.app._t("dlg.files.title"),
                                       self.app._t("dlg.files.nothing"))
                return
            self.app.stop_work_panel()
            self._ok(selected["files"])

        self.app.work_panel_footer(panel, on_ok=finish,
                                   ok_text=self.app._t("btn.continue"),
                                   cancel=True,
                                   cancel_text=self.app._t("btn.cancel"))

    def _dialog_channel_selector(self, kw):
        mips = kw["mips"]
        n = len(mips)
        if n < 1:
            self._ok({
                "selected_channels": [0],
                "cell_idx": 0,
                "has_dapi": False,
                "dapi_idx": None,
            })
            return

        cell_default = kw.get("cell_idx_default", 0)
        dapi_default = kw.get("dapi_idx_default", None)
        has_dapi_default = kw.get("has_dapi_default", False)

        if not (0 <= int(cell_default) < n):
            cell_default = 0
        if dapi_default is not None and not (0 <= int(dapi_default) < n):
            dapi_default = None
        if dapi_default is not None and int(dapi_default) == int(cell_default):
            dapi_default = None
            has_dapi_default = False

        panel = self.app.start_work_panel(self.app._t("dlg.channels.title"))

        tk.Label(panel,
                 text=self.app._t("dlg.channels.instructions"),
                 justify="left", anchor="w",
                 bg=self.app.BG_CARD, fg=self.app.FG,
                 font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(0, 10))

        var_no_dapi = tk.BooleanVar(value=(not has_dapi_default))
        chk_no_dapi = tk.Checkbutton(
            panel, text=self.app._t("dlg.channels.no_dapi_checkbox"),
            variable=var_no_dapi,
            bg=self.app.BG_CARD, fg=self.app.FG,
            activebackground=self.app.BG_CARD,
            activeforeground=self.app.ACCENT,
            selectcolor=self.app.BG_SOFT,
            relief="flat", bd=0,
            font=("Segoe UI", 10), cursor="hand2",
            highlightthickness=0,
        )
        chk_no_dapi.pack(anchor="w", pady=(0, 10))

        scroll_frame = tk.Frame(panel, bg=self.app.BG_CARD)
        scroll_frame.pack(fill="both", expand=True)

        canvas = tk.Canvas(scroll_frame, bg=self.app.BG_CARD, highlightthickness=0)

        scrollbar = RoundedScrollbar(scroll_frame, command=canvas.xview,
                                     width=10, orient="horizontal", bg=self.app.BG_CARD,
                                     thumb=self.app.BORDER, thumb_hover=self.app.ACCENT)

        canvas.configure(xscrollcommand=scrollbar.set)

        scrollbar.pack(side="bottom", fill="x", pady=(4, 0))
        canvas.pack(side="top", fill="both", expand=True)

        imgs_frame = tk.Frame(canvas, bg=self.app.BG_CARD)
        canvas.create_window((0, 0), window=imgs_frame, anchor="nw")

        def _on_frame_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        imgs_frame.bind("<Configure>", _on_frame_configure)

        def _on_mousewheel(event):
            canvas.xview_scroll(int(-1 * (event.delta / 120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        selected_flags = []
        for i in range(n):
            if has_dapi_default:
                selected_flags.append(i in (int(cell_default), int(dapi_default)))
            else:
                selected_flags.append(i == int(cell_default))

        selected_vars = [tk.BooleanVar(value=selected_flags[i]) for i in range(n)]
        cell_choice = {"val": int(cell_default)}
        dapi_choice = {"val": int(dapi_default) if has_dapi_default else None}
        tile_drawers = {}
        tile_photos = {}

        def _refresh_all():
            for f in tile_drawers.values():
                f()

        def _pick_cell(index):
            cell_choice["val"] = index
            if not var_no_dapi.get() and dapi_choice["val"] == index:
                for k in range(n):
                    if k != index:
                        dapi_choice["val"] = k
                        break
            _refresh_all()

        def _pick_dapi(index):
            if var_no_dapi.get():
                return
            dapi_choice["val"] = index
            if cell_choice["val"] == index:
                for k in range(n):
                    if k != index:
                        cell_choice["val"] = k
                        break
            _refresh_all()

        def _toggle_selected(index):
            selected_vars[index].set(not selected_vars[index].get())
            if not selected_vars[index].get():
                if cell_choice["val"] == index:
                    for k in range(n):
                        if k != index and selected_vars[k].get():
                            cell_choice["val"] = k
                            break
                if dapi_choice["val"] == index:
                    for k in range(n):
                        if k != index and selected_vars[k].get():
                            dapi_choice["val"] = k
                            break
            _refresh_all()

        def _make_tile(index, mip):
            arr = self._to_uint8(mip)
            pil = self._pil_from_gray(arr)
            w, h, _ = self._fit_size(pil.width, pil.height,
                                     max_w=200, max_h=200)
            pil = pil.resize((w, h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)
            tile_photos[index] = photo

            tile_w = w + 8
            tile_h = h + 8 + 90

            tile = tk.Canvas(imgs_frame, width=tile_w, height=tile_h,
                             bg=self.app.BG_CARD,
                             highlightthickness=0, bd=0)
            tile.pack(side="left", padx=8)
            tile.image = photo

            def draw_tile():
                tile.delete("all")
                is_selected = selected_vars[index].get()
                is_cell = (cell_choice["val"] == index)
                is_dapi = (dapi_choice["val"] == index) and (not var_no_dapi.get())

                if is_cell and is_selected:
                    border = self.app.OK
                elif is_dapi and is_selected:
                    border = self.app.ACCENT
                elif is_selected:
                    border = self.app.BORDER
                else:
                    border = self.app.FG_MUTED

                tile.create_rectangle(1, 1, tile_w - 1, tile_h - 1,
                                      outline=border, width=2,
                                      fill=self.app.BG_SOFT)
                tile.create_image(4, 4, anchor="nw", image=photo)

                tile.create_text(tile_w // 2, tile_h - 72,
                                 text=self.app._t("dlg.channels.channel_n", n=index),
                                 fill=self.app.FG,
                                 font=("Segoe UI", 10, "bold"))

                if is_cell:
                    role_lbl = self.app._t("dlg.channels.role_cells")
                    role_col = self.app.OK
                elif is_dapi:
                    role_lbl = self.app._t("dlg.channels.role_dapi")
                    role_col = self.app.ACCENT
                else:
                    role_lbl = self.app._t("dlg.channels.role_none")
                    role_col = self.app.FG_DIM
                tile.create_text(tile_w // 2, tile_h - 54,
                                 text=role_lbl,
                                 fill=role_col,
                                 font=("Segoe UI", 9, "bold"))

                cb_x1 = tile_w // 2 - 60
                cb_x2 = tile_w // 2 + 60
                cb_y1 = tile_h - 36
                cb_y2 = tile_h - 12
                tile.create_rectangle(cb_x1, cb_y1, cb_x2, cb_y2,
                                      outline=self.app.BORDER, width=1,
                                      fill=self.app.BG_CARD)
                tile.create_text((cb_x1 + cb_x2) // 2, (cb_y1 + cb_y2) // 2,
                                 text=(self.app._t("dlg.channels.take")
                                       if is_selected
                                       else self.app._t("dlg.channels.skip")),
                                 fill=(self.app.OK if is_selected
                                       else self.app.FG_DIM),
                                 font=("Segoe UI", 9))

            def on_left(event, idx=index, th=tile_h):
                if event.y > th - 42:
                    _toggle_selected(idx)
                else:
                    _pick_cell(idx)

            def on_right(event, idx=index):
                _pick_dapi(idx)

            tile.bind("<Button-1>", on_left)
            tile.bind("<Button-3>", on_right)

            tile_drawers[index] = draw_tile
            draw_tile()

        for i in range(n):
            _make_tile(i, mips[i])

        def _on_no_dapi_toggle():
            if var_no_dapi.get():
                dapi_choice["val"] = None
            else:
                for k in range(n):
                    if k != cell_choice["val"] and selected_vars[k].get():
                        dapi_choice["val"] = k
                        break
            _refresh_all()

        chk_no_dapi.configure(command=_on_no_dapi_toggle)

        def finish():
            canvas.unbind_all("<MouseWheel>")
            selected = [i for i in range(n) if selected_vars[i].get()]
            if not selected:
                messagebox.showwarning(self.app._t("dlg.channels.title"),
                                       self.app._t("dlg.channels.choose_one"))
                return
            if cell_choice["val"] not in selected:
                selected_vars[cell_choice["val"]].set(True)
                selected = [i for i in range(n) if selected_vars[i].get()]

            has_dapi = not var_no_dapi.get()
            dapi_idx = dapi_choice["val"] if has_dapi else None

            if has_dapi and (dapi_idx is None or dapi_idx not in selected):
                alt = [i for i in selected if i != cell_choice["val"]]
                if alt:
                    dapi_idx = alt[0]
                else:
                    has_dapi = False
                    dapi_idx = None

            self.app.stop_work_panel()
            self._ok({
                "selected_channels": selected,
                "cell_idx": int(cell_choice["val"]),
                "has_dapi": bool(has_dapi),
                "dapi_idx": int(dapi_idx) if dapi_idx is not None else None,
            })

        def cancel():
            canvas.unbind_all("<MouseWheel>")
            self.app.stop_work_panel()
            self._ok(None)

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            cancel=True, cancel_value=None,
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.channels.help")
        )

        for child in self.app._footer_btns.winfo_children():
            if isinstance(child, RoundedButton) and child._label == self.app._t("btn.cancel"):
                child._command = cancel

    def _dialog_dapi_channel(self, kw):
        mips = kw["mips"]
        if len(mips) < 2:
            self._ok(0); return

        panel = self.app.start_work_panel(self.app._t("dlg.dapi.title"))

        tk.Label(panel,
                 text=self.app._t("dlg.dapi.text"),
                 justify="left", anchor="w",
                 bg=self.app.BG_CARD, fg=self.app.FG,
                 font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(0, 10))

        imgs_frame = tk.Frame(panel, bg=self.app.BG_CARD)
        imgs_frame.pack(anchor="w")

        result = {"val": None}

        def make_tile(index, mip):
            arr = self._to_uint8(mip)
            pil = self._pil_from_gray(arr)
            w, h, _ = self._fit_size(pil.width, pil.height, max_w=240, max_h=260)
            pil = pil.resize((w, h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)

            tile_w = w + 8
            tile_h = h + 8 + 28

            tile = tk.Canvas(imgs_frame, width=tile_w, height=tile_h,
                             bg=self.app.BG_CARD,
                             highlightthickness=0, bd=0)
            tile.pack(side="left", padx=8)
            tile.image = photo

            def draw_tile(hover=False):
                tile.delete("all")
                r = 12
                border = self.app.ACCENT if hover else self.app.BORDER
                fill = self.app.BG_SOFT
                pts = [
                    1 + r, 1, tile_w - 1 - r, 1, tile_w - 1, 1,
                    tile_w - 1, 1 + r, tile_w - 1, tile_h - 1 - r,
                    tile_w - 1, tile_h - 1,
                    tile_w - 1 - r, tile_h - 1, 1 + r, tile_h - 1,
                    1, tile_h - 1, 1, tile_h - 1 - r,
                    1, 1 + r, 1, 1,
                    ]
                tile.create_polygon(pts, smooth=True,
                                    fill=fill, outline=border, width=2)
                tile.create_image(4, 4, anchor="nw", image=photo)
                tile.create_text(tile_w // 2, tile_h - 12,
                                 text=self.app._t("dlg.channels.channel_n", n=index),
                                 fill=self.app.FG,
                                 font=("Segoe UI", 10, "bold"))

            draw_tile(False)

            def on_enter(_):
                draw_tile(True)
                tile.configure(cursor="hand2")

            def on_leave(_):
                draw_tile(False)
                tile.configure(cursor="")

            tile.bind("<Enter>", on_enter)
            tile.bind("<Leave>", on_leave)
            tile.bind("<Button-1>", lambda e, i=index: on_pick(i))

        def on_pick(i):
            result["val"] = i
            self.app.stop_work_panel()
            self._ok(i)

        for i in range(2):
            make_tile(i, mips[i])

        self.app.work_panel_footer(panel, on_ok=None, ok_text=None,
                                   cancel=True, cancel_value=None,
                                   cancel_text=self.app._t("btn.cancel"),
                                   help_text=self.app._t("dlg.dapi.help"))

    def _build_color_mip(self, astro_norm):
        z_dim, y_dim, x_dim = astro_norm.shape

        def _norm_channel(ch):
            nz = ch[ch > 0]
            if nz.size == 0:
                return np.zeros_like(ch, dtype=np.float32)
            p0, p99 = np.percentile(nz, [1, 99])
            if p99 - p0 < 1e-6:
                return np.zeros_like(ch, dtype=np.float32)
            return np.clip((ch - p0) / (p99 - p0), 0.0, 1.0).astype(np.float32)

        if z_dim < 3:
            mip = np.max(astro_norm, axis=0)
            mip_n = _norm_channel(mip)
            rgb = np.stack([mip_n, mip_n, mip_n], axis=-1)
            return (rgb * 255).astype(np.uint8)

        third = z_dim // 3

        r = _norm_channel(np.max(astro_norm[0:third, :, :], axis=0))
        g = _norm_channel(np.max(astro_norm[third:2 * third, :, :], axis=0))
        b = _norm_channel(np.max(astro_norm[2 * third:, :, :], axis=0))

        mip_gray = np.maximum(np.maximum(r, g), b)

        alpha = 0.3

        r_out = mip_gray * (1.0 - alpha) + r * alpha
        g_out = mip_gray * (1.0 - alpha) + g * alpha
        b_out = mip_gray * (1.0 - alpha) + b * alpha

        rgb = np.stack([r_out, g_out, b_out], axis=-1)
        rgb = np.clip(rgb, 0.0, 1.0)
        return (rgb * 255).astype(np.uint8)

    def _dialog_soma_confirm(self, kw):
        astro_norm = kw["astro_norm"]
        click_centers = kw["click_centers"]

        panel = self.app.start_work_panel(self.app._t("dlg.soma_confirm.title"))

        z_dim, img_h, img_w = astro_norm.shape
        z_var = tk.IntVar(value=0)
        photo_ref = {"photo": None}

        slider_holder = tk.Frame(panel, bg=self.app.BG_CARD)
        slider_holder.pack(side="bottom", fill="x", pady=(8, 0))

        head = tk.Frame(slider_holder, bg=self.app.BG_CARD)
        head.pack(fill="x")
        tk.Label(head, text=self.app._t("dlg.soma_confirm.z_slice"),
                 bg=self.app.BG_CARD, fg=self.app.FG_DIM,
                 font=("Segoe UI", 9, "bold")).pack(side="left")
        z_value_lbl = tk.Label(head, text="0",
                               bg=self.app.BG_CARD, fg=self.app.FG,
                               font=("Consolas", 10, "bold"))
        z_value_lbl.pack(side="right")

        info = tk.Label(panel, text="",
                        bg=self.app.BG_CARD, fg=self.app.FG_DIM,
                        font=("Segoe UI", 9), anchor="w")
        info.pack(side="bottom", anchor="w", pady=(4, 0))

        canvas = tk.Canvas(panel, bg="#101216", highlightthickness=0)
        canvas.pack(side="top", fill="both", expand=True, pady=(4, 8))

        def redraw_canvas(event=None):
            z = z_var.get()
            cw = canvas.winfo_width()
            ch = canvas.winfo_height()
            if cw <= 2 or ch <= 2:
                return
            s = min(cw / img_w, ch / img_h)
            new_w = max(1, int(img_w * s))
            new_h = max(1, int(img_h * s))

            pil = self._pil_from_gray(self._to_uint8(astro_norm[z])).resize(
                (new_w, new_h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)
            photo_ref["photo"] = photo

            canvas.delete("all")
            x0 = (cw - new_w) // 2
            y0 = (ch - new_h) // 2
            canvas.create_image(x0, y0, anchor="nw", image=photo)

            scale_x = new_w / img_w
            scale_y = new_h / img_h
            n_here = 0
            for (cz, cy, cx) in click_centers:
                if cz == z:
                    n_here += 1
                    px = x0 + cx * scale_x
                    py = y0 + cy * scale_y
                    r_ = 6
                    canvas.create_oval(px - r_, py - r_, px + r_, py + r_,
                                       outline="white", fill="#8fd694", width=2)

            badge_w, badge_h = 260, 30
            bx, by = x0 + 8, y0 + 8
            canvas.create_rectangle(bx, by, bx + badge_w, by + badge_h,
                                    outline="#8fd694", width=2, fill="#101216")
            canvas.create_text(bx + 10, by + badge_h // 2, anchor="w",
                               text=self.app._t("dlg.soma_confirm.z_badge",
                                                z=z, n=n_here),
                               fill="#8fd694", font=("Segoe UI", 11, "bold"))

            info.config(text=self.app._t("dlg.soma_confirm.info",
                                         total=len(click_centers),
                                         here=n_here))

        canvas.bind("<Configure>", lambda e: redraw_canvas())

        def on_slide(v):
            z_value_lbl.config(text=str(int(float(v))))
            redraw_canvas()

        RoundedSlider(slider_holder,
                      from_=0, to=z_dim - 1,
                      variable=z_var,
                      length=200,
                      command=on_slide,
                      bg=self.app.BG_CARD,
                      track=self.app.BORDER,
                      fill=self.app.ACCENT,
                      thumb="#d4d4dc",
                      thumb_hover="#ffffff",
                      fg=self.app.FG,
                      showvalue=False).pack(fill="x", pady=(2, 0))

        def finish():
            self.app.stop_work_panel()
            self._ok(True)

        canvas.after(30, redraw_canvas)

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.start_analysis"),
            cancel=True, cancel_value=False,
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.soma_confirm.help")
        )

    def _dialog_soma_picker(self, kw):
        mip = kw["mip"]
        astro_norm = kw.get("astro_norm", None)

        panel = self.app.start_work_panel(self.app._t("dlg.soma.title"))

        if astro_norm is not None:
            color_mip = self._build_color_mip(astro_norm)
        else:
            u8 = self._to_uint8(mip)
            color_mip = np.stack([u8, u8, u8], axis=-1)

        img_h, img_w = color_mip.shape[:2]

        points = []
        scale_ref = {"s": 1.0, "x0": 0, "y0": 0}
        photo_ref = {"photo": None}
        cursor_ref = {"x": None, "y": None}

        OV_IMG = 260
        OV_PAD = 8

        overlay_ref = {
            "active": False,
            "frame": None,
            "canvas": None,
            "photo": None,
            "z_var": None,
            "z_lbl": None,
            "crop": None,
            "crop_origin": None,
            "crop_size": None,
        }

        counter = tk.Label(panel, text=self.app._t("dlg.soma.counter", n=0),
                           bg=self.app.BG_CARD, fg=self.app.FG_DIM,
                           font=("Segoe UI", 9))
        counter.pack(side="bottom", anchor="w", pady=(6, 0))

        canvas = tk.Canvas(panel, bg="#101216", highlightthickness=0)
        canvas.pack(side="top", fill="both", expand=True, pady=4)

        def _draw_crosshair():
            canvas.delete("crosshair")
            if overlay_ref["active"]:
                return
            x = cursor_ref["x"]
            y = cursor_ref["y"]
            if x is None or y is None:
                return
            cw = canvas.winfo_width()
            ch = canvas.winfo_height()
            if cw <= 2 or ch <= 2:
                return

            canvas.create_line(0, y, cw, y,
                               fill="#ffd54a", width=1, tags="crosshair")
            canvas.create_line(x, 0, x, ch,
                               fill="#ffd54a", width=1, tags="crosshair")
            r = 3
            canvas.create_oval(x - r, y - r, x + r, y + r,
                               outline="#1a1c22", fill="#ffd54a", width=1,
                               tags="crosshair")

            coords_text = _get_cursor_orig_text(x, y)
            if coords_text is not None:
                canvas.create_text(10, 10, anchor="nw",
                                   text=coords_text,
                                   fill="#ffd54a",
                                   font=("Consolas", 10, "bold"),
                                   tags="crosshair")

        def _get_cursor_orig_text(sx, sy):
            s = scale_ref["s"]
            x0 = scale_ref["x0"]
            y0 = scale_ref["y0"]
            if s <= 0:
                return None
            ox = int((sx - x0) / s)
            oy = int((sy - y0) / s)
            if not (0 <= ox < img_w and 0 <= oy < img_h):
                return None
            return self.app._t("dlg.soma.coords", y=oy, x=ox)

        def _draw_overlay():
            if not overlay_ref["active"]:
                return
            ov_canvas = overlay_ref["canvas"]
            if ov_canvas is None:
                return

            crop = overlay_ref["crop"]
            z_var = overlay_ref["z_var"]
            z_lbl = overlay_ref["z_lbl"]
            crop_origin = overlay_ref["crop_origin"]
            crop_size = overlay_ref["crop_size"]
            if (crop is None or z_var is None or crop_origin is None
                    or crop_size is None):
                return

            z_cur = int(z_var.get())
            z_dim_crop = crop.shape[0]
            if z_cur < 0:
                z_cur = 0
            if z_cur >= z_dim_crop:
                z_cur = z_dim_crop - 1

            sl = crop[z_cur]
            nz = sl[sl > 0]
            if nz.size > 0:
                p2, p98 = np.percentile(nz, [2, 98])
                if p98 - p2 < 1e-6:
                    sl_n = np.zeros_like(sl)
                else:
                    sl_n = np.clip((sl - p2) / (p98 - p2), 0.0, 1.0)
            else:
                sl_n = np.zeros_like(sl)

            u8 = (sl_n * 255).astype(np.uint8)
            pil = self._pil_from_gray(u8)
            pil_resized = pil.resize((OV_IMG, OV_IMG), Image.BILINEAR)

            base_rgb = np.array(pil_resized.convert("RGB"))

            origin_y, origin_x = crop_origin
            half = crop_size // 2

            for pt in points:
                py, px = pt[0], pt[1]
                dy = py - origin_y + half
                dx = px - origin_x + half
                if 0 <= dy < crop_size and 0 <= dx < crop_size:
                    cy_ = int(dy * OV_IMG / crop_size)
                    cx_ = int(dx * OV_IMG / crop_size)
                    r_ = 4
                    y0_ = max(0, cy_ - r_)
                    y1_ = min(OV_IMG, cy_ + r_)
                    x0_ = max(0, cx_ - r_)
                    x1_ = min(OV_IMG, cx_ + r_)
                    if y1_ > y0_ and x1_ > x0_:
                        base_rgb[y0_:y1_, x0_:x1_, 0] = 143
                        base_rgb[y0_:y1_, x0_:x1_, 1] = 214
                        base_rgb[y0_:y1_, x0_:x1_, 2] = 148

            photo = ImageTk.PhotoImage(Image.fromarray(base_rgb, mode="RGB"))
            overlay_ref["photo"] = photo

            ov_canvas.delete("all")
            ov_canvas.create_image(0, 0, anchor="nw", image=photo)

            if z_lbl is not None:
                z_lbl.config(text=self.app._t("dlg.soma.overlay_z", z=z_cur))

        def _close_overlay():
            if not overlay_ref["active"]:
                return
            if overlay_ref["frame"] is not None:
                overlay_ref["frame"].destroy()
                overlay_ref["frame"] = None
            overlay_ref["canvas"] = None
            overlay_ref["photo"] = None
            overlay_ref["z_var"] = None
            overlay_ref["z_lbl"] = None
            overlay_ref["crop"] = None
            overlay_ref["crop_origin"] = None
            overlay_ref["crop_size"] = None
            overlay_ref["active"] = False
            redraw()

        def _on_overlay_slide(_v=None):
            _draw_overlay()

        def _on_overlay_click(event):
            if not overlay_ref["active"]:
                return
            z_var = overlay_ref["z_var"]
            crop_origin = overlay_ref["crop_origin"]
            crop_size = overlay_ref["crop_size"]
            if (z_var is None or crop_origin is None or crop_size is None):
                return

            if not (0 <= event.x < OV_IMG and 0 <= event.y < OV_IMG):
                return

            cx_in_crop = int(event.x * crop_size / OV_IMG)
            cy_in_crop = int(event.y * crop_size / OV_IMG)

            origin_y, origin_x = crop_origin
            half = crop_size // 2

            orig_y = origin_y - half + cy_in_crop
            orig_x = origin_x - half + cx_in_crop

            if not (0 <= orig_y < img_h and 0 <= orig_x < img_w):
                return

            z_cur = int(z_var.get())

            replaced = False
            for idx, pt in enumerate(points):
                if pt[0] == orig_y and pt[1] == orig_x:
                    points[idx] = (orig_y, orig_x, z_cur)
                    replaced = True
                    break
            if not replaced:
                points.append((orig_y, orig_x, z_cur))

            counter.config(text=self.app._t("dlg.soma.counter", n=len(points)))
            redraw()
            _draw_overlay()

        def redraw():
            cw = canvas.winfo_width()
            ch = canvas.winfo_height()
            if cw <= 2 or ch <= 2:
                return

            s = min(cw / img_w, ch / img_h)
            new_w = max(1, int(img_w * s))
            new_h = max(1, int(img_h * s))
            x0 = (cw - new_w) // 2
            y0 = (ch - new_h) // 2

            scale_ref["s"] = s
            scale_ref["x0"] = x0
            scale_ref["y0"] = y0

            pil = Image.fromarray(color_mip, mode="RGB").resize(
                (new_w, new_h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)
            photo_ref["photo"] = photo

            canvas.delete("all")
            canvas.create_image(x0, y0, anchor="nw", image=photo)

            r = 6
            for pt in points:
                py, px = pt[0], pt[1]
                dx = x0 + px * s
                dy = y0 + py * s
                canvas.create_oval(dx - r, dy - r, dx + r, dy + r,
                                   outline="white", fill="#8fd694", width=2)

            _draw_crosshair()

        canvas.bind("<Configure>", lambda e: redraw())

        def _get_orig_coords(event):
            s = scale_ref["s"]
            x0 = scale_ref["x0"]
            y0 = scale_ref["y0"]
            if s <= 0:
                return None
            x = int((event.x - x0) / s)
            y = int((event.y - y0) / s)
            if not (0 <= x < img_w and 0 <= y < img_h):
                return None
            return y, x

        def _auto_pick_z(orig_y, orig_x):
            z_cur = 0
            best_score = -1.0
            if astro_norm is not None:
                for zz in range(astro_norm.shape[0]):
                    win = astro_norm[zz,
                    max(0, orig_y - 2):min(img_h, orig_y + 3),
                    max(0, orig_x - 2):min(img_w, orig_x + 3)]
                    if win.size > 0:
                        sc = float(np.mean(win))
                        if sc > best_score:
                            best_score = sc
                            z_cur = zz
            return z_cur

        def on_motion(event):
            cursor_ref["x"] = event.x
            cursor_ref["y"] = event.y
            _draw_crosshair()

        def on_leave(event):
            cursor_ref["x"] = None
            cursor_ref["y"] = None
            canvas.delete("crosshair")

        def on_left(event):
            if overlay_ref["active"]:
                return
            coords = _get_orig_coords(event)
            if coords is None:
                return
            y, x = coords

            ctrl_held = bool(event.state & 0x0004)
            if ctrl_held:
                if points:
                    d = [(pt[0] - y) ** 2 + (pt[1] - x) ** 2 for pt in points]
                    points.pop(int(np.argmin(d)))
            else:
                z_cur = _auto_pick_z(y, x)
                points.append((y, x, z_cur))

            counter.config(text=self.app._t("dlg.soma.counter", n=len(points)))
            redraw()

        def on_right(event):
            if overlay_ref["active"]:
                return
            coords = _get_orig_coords(event)
            if coords is None:
                return
            y, x = coords
            _open_overlay(y, x)

        canvas.bind("<Button-1>", on_left)
        canvas.bind("<Button-3>", on_right)
        canvas.bind("<Motion>", on_motion)
        canvas.bind("<Leave>", on_leave)

        def _open_overlay(orig_y, orig_x):
            _close_overlay()

            if astro_norm is None:
                return

            w = max(40, min(img_h, img_w) // 10)
            side = 2 * w

            crop = np.zeros((astro_norm.shape[0], side, side),
                            dtype=astro_norm.dtype)

            y0c = orig_y - w
            y1c = orig_y + w
            x0c = orig_x - w
            x1c = orig_x + w

            sy0 = max(0, -y0c)
            sy1 = sy0 + (min(y1c, img_h) - max(y0c, 0))
            sx0 = max(0, -x0c)
            sx1 = sx0 + (min(x1c, img_w) - max(x0c, 0))

            dy0 = max(0, y0c)
            dy1 = dy0 + (sy1 - sy0)
            dx0 = max(0, x0c)
            dx1 = dx0 + (sx1 - sx0)

            crop[:, sy0:sy1, sx0:sx1] = astro_norm[:, dy0:dy1, dx0:dx1]

            crop_size = side
            crop_origin = (orig_y, orig_x)

            ov_frame = tk.Frame(canvas, bg=self.app.BG_CARD,
                                highlightthickness=2,
                                highlightbackground=self.app.ACCENT)

            ov_canvas = tk.Canvas(ov_frame, bg="#101216",
                                  highlightthickness=0, bd=0,
                                  width=OV_IMG, height=OV_IMG)
            ov_canvas.pack(side="top", padx=0, pady=0)

            z_max = crop.shape[0] - 1
            z_best = 0
            best_score = -1.0
            for zz in range(crop.shape[0]):
                sl = crop[zz]
                if sl.size == 0:
                    continue
                score = float(np.max(sl))
                if score > best_score:
                    best_score = score
                    z_best = zz

            z_var = tk.IntVar(value=int(z_best))

            slider_holder = tk.Frame(ov_frame, bg=self.app.BG_CARD)
            slider_holder.pack(side="top", fill="x", padx=OV_PAD, pady=(4, 0))

            z_lbl = tk.Label(slider_holder,
                             text=self.app._t("dlg.soma.overlay_z", z=int(z_best)),
                             bg=self.app.BG_CARD, fg=self.app.FG,
                             font=("Consolas", 9, "bold"))
            z_lbl.pack(side="top", anchor="w")

            RoundedSlider(slider_holder,
                          from_=0, to=max(0, z_max),
                          variable=z_var,
                          length=OV_IMG - 20,
                          command=_on_overlay_slide,
                          bg=self.app.BG_CARD,
                          track=self.app.BORDER,
                          fill=self.app.ACCENT,
                          thumb="#d4d4dc",
                          thumb_hover="#ffffff",
                          fg=self.app.FG,
                          showvalue=False).pack(fill="x", pady=(2, 0))

            btn_holder = tk.Frame(ov_frame, bg=self.app.BG_CARD)
            btn_holder.pack(side="bottom", fill="x", padx=OV_PAD, pady=(4, 6))

            RoundedButton(btn_holder, self.app._t("btn.ok"),
                          command=_close_overlay,
                          width=80, height=26, radius=13,
                          bg=self.app.ACCENT, fg="#1a1c22",
                          hover_bg=self.app.ACCENT_D, hover_fg="#ffffff",
                          active_bg="#4a6fbf",
                          font=("Segoe UI", 9, "bold")).pack(side="right")

            ov_frame.update_idletasks()
            need_w = ov_frame.winfo_reqwidth()
            need_h = ov_frame.winfo_reqheight()

            cw = canvas.winfo_width()

            ov_x = cw - need_w - 10
            ov_y = 10
            if ov_x < 0:
                ov_x = 0

            ov_frame.place(x=ov_x, y=ov_y, width=need_w, height=need_h)
            ov_frame.update_idletasks()

            overlay_ref["active"] = True
            overlay_ref["frame"] = ov_frame
            overlay_ref["canvas"] = ov_canvas
            overlay_ref["z_var"] = z_var
            overlay_ref["z_lbl"] = z_lbl
            overlay_ref["crop"] = crop
            overlay_ref["crop_origin"] = crop_origin
            overlay_ref["crop_size"] = crop_size

            ov_canvas.bind("<Button-1>", _on_overlay_click)

            _draw_overlay()
            _draw_crosshair()

        def finish():
            _close_overlay()
            self.app.stop_work_panel()
            self._ok(list(points))

        canvas.after(30, redraw)

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.done"), cancel=True,
            cancel_value=[],
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.soma.help")
        )

    def _dialog_slice_selector(self, kw):
        stack = kw["stack"]
        z_dim = stack.shape[0]
        title = kw.get("title") or self.app._t("dlg.slices.title")
        panel = self.app.start_work_panel(title)

        u8_stack = np.stack([self._to_uint8(stack[z]) for z in range(z_dim)], axis=0)
        img_h, img_w = u8_stack.shape[1:]

        status = {z: True for z in range(z_dim)}
        z_var = tk.IntVar(value=0)
        photo_ref = {"img": None, "photo": None, "scale": 1.0}

        slider_holder = tk.Frame(panel, bg=self.app.BG_CARD)
        slider_holder.pack(side="bottom", fill="x", pady=(8, 0))

        head = tk.Frame(slider_holder, bg=self.app.BG_CARD)
        head.pack(fill="x")
        tk.Label(head, text=self.app._t("dlg.slices.z_slice"),
                 bg=self.app.BG_CARD, fg=self.app.FG_DIM,
                 font=("Segoe UI", 9, "bold")).pack(side="left")
        z_value_lbl = tk.Label(head, text="0",
                               bg=self.app.BG_CARD, fg=self.app.FG,
                               font=("Consolas", 10, "bold"))
        z_value_lbl.pack(side="right")

        toggle_bar = tk.Frame(panel, bg=self.app.BG_CARD)
        toggle_bar.pack(side="bottom", fill="x", pady=(6, 0))

        canvas = tk.Canvas(panel, bg="#101216", highlightthickness=0)
        canvas.pack(side="top", fill="both", expand=True, pady=(4, 8))

        def redraw_canvas(event=None):
            z = z_var.get()
            cw = canvas.winfo_width()
            ch = canvas.winfo_height()
            if cw <= 2 or ch <= 2:
                return

            s = min(cw / img_w, ch / img_h)
            new_w = max(1, int(img_w * s))
            new_h = max(1, int(img_h * s))

            photo_ref["scale"] = s
            pil = self._pil_from_gray(u8_stack[z]).resize(
                (new_w, new_h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)
            photo_ref["photo"] = photo

            canvas.delete("all")
            x0 = (cw - new_w) // 2
            y0 = (ch - new_h) // 2
            canvas.create_image(x0, y0, anchor="nw", image=photo)

            color = "#8fd694" if status[z] else "#e07a7a"
            label = (self.app._t("dlg.slices.kept") if status[z]
                     else self.app._t("dlg.slices.deleted"))
            badge_w, badge_h = 220, 30
            bx = x0 + 8
            by = y0 + 8
            canvas.create_rectangle(bx, by, bx + badge_w, by + badge_h,
                                    outline=color, width=2, fill="#101216")
            canvas.create_text(bx + 10, by + badge_h // 2, anchor="w",
                               text=self.app._t("dlg.slices.badge",
                                                z=z, label=label),
                               fill=color, font=("Segoe UI", 11, "bold"))

        canvas.bind("<Configure>", lambda e: redraw_canvas())

        def render():
            z_var.set(z_var.get())
            z_value_lbl.config(text=str(z_var.get()))
            redraw_canvas()

        def toggle():
            z = z_var.get()
            status[z] = not status[z]
            render()

        RoundedButton(toggle_bar, self.app._t("dlg.slices.toggle"),
                      command=toggle,
                      width=220, height=32, radius=16,
                      bg=self.app.BG_SOFT, fg=self.app.FG,
                      hover_bg=self.app.ACCENT, hover_fg="#1a1c22",
                      active_bg=self.app.ACCENT_D,
                      font=("Segoe UI", 9)).pack(anchor="w")

        RoundedSlider(slider_holder,
                      from_=0, to=z_dim - 1,
                      variable=z_var,
                      length=200,
                      command=lambda *_: render(),
                      bg=self.app.BG_CARD,
                      track=self.app.BORDER,
                      fill=self.app.ACCENT,
                      thumb="#d4d4dc",
                      thumb_hover="#ffffff",
                      fg=self.app.FG,
                      showvalue=False).pack(fill="x", pady=(2, 0))

        def finish():
            self.app.stop_work_panel()
            self._ok([z for z in range(z_dim) if not status[z]])

        canvas.after(30, redraw_canvas)
        z_value_lbl.config(text="0")

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            cancel=True, cancel_value=[],
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.slices.help")
        )

    def _dialog_cell_filter(self, kw):
        mip, cell_mip, colors = kw["mip"], kw["cell_mip"], kw["colors"]
        y_dim, x_dim = mip.shape
        panel = self.app.start_work_panel(self.app._t("dlg.filter.title"))

        base_u8 = self._to_uint8(mip)
        disp_w, disp_h, scale = self._fit_size(x_dim, y_dim, max_w=500, max_h=470)
        base_small = np.array(
            self._pil_from_gray(base_u8).resize((disp_w, disp_h), Image.BILINEAR))

        canvas = tk.Canvas(panel, width=disp_w, height=disp_h, bg="black",
                           highlightthickness=0)
        canvas.pack(pady=4)

        state = {int(c): 0 for c in np.unique(cell_mip)[1:]}
        photo_ref = {"img": None}

        masks_small = {}
        for cid in state:
            m = (cell_mip == cid)
            if not np.any(m):
                continue
            small = np.array(
                Image.fromarray(m.astype(np.uint8) * 255, mode="L")
                .resize((disp_w, disp_h), Image.NEAREST)) > 127
            masks_small[cid] = small

        selected_colors = {1: (255, 80, 80), 2: (80, 130, 255)}

        def render():
            rgb = np.stack([base_small, base_small, base_small], axis=-1).astype(np.float32)
            for cid, small in masks_small.items():
                if not small.any():
                    continue
                st = state[cid]
                if st == 0:
                    col = np.array(colors.get(cid, (1.0, 0.5, 0.0)),
                                   dtype=np.float32) * 255.0
                    a = 0.45
                else:
                    col = np.array(selected_colors[st], dtype=np.float32)
                    a = 0.7
                rgb[small] = rgb[small] * (1 - a) + col * a
            pil = Image.fromarray(rgb.astype(np.uint8), mode="RGB")
            photo = ImageTk.PhotoImage(pil)
            photo_ref["img"] = photo
            canvas.delete("all")
            canvas.create_image(0, 0, anchor="nw", image=photo)

        def on_click(event):
            x = int(event.x / scale)
            y = int(event.y / scale)
            if not (0 <= x < x_dim and 0 <= y < y_dim):
                return
            cid = int(cell_mip[y, x])
            if cid not in state:
                return
            if state[cid] != 0:
                state[cid] = 0
            elif event.num == 1:
                state[cid] = 2
            elif event.num == 3:
                state[cid] = 1
            render()

        canvas.bind("<Button-1>", on_click)
        canvas.bind("<Button-3>", on_click)
        render()

        def finish():
            self.app.stop_work_panel()
            self._ok({"skip": False,
                      "states": {c: s for c, s in state.items() if s != 0}})

        def skip():
            self.app.stop_work_panel()
            self._ok({"skip": True, "states": {}})

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.continue"),
            extra=[(self.app._t("btn.skip"), skip)],
            help_text=self.app._t("dlg.filter.help")
        )

    def _dialog_resegment_review(self, kw):
        mip = kw["mip"]; new_masks = kw["new_masks"]
        y_dim, x_dim = mip.shape
        panel = self.app.start_work_panel(self.app._t("dlg.review.title"))

        base_u8 = self._to_uint8(mip)
        disp_w, disp_h, scale = self._fit_size(x_dim, y_dim, max_w=500, max_h=470)
        base_small = np.array(
            self._pil_from_gray(base_u8).resize((disp_w, disp_h), Image.BILINEAR))

        canvas = tk.Canvas(panel, width=disp_w, height=disp_h, bg="black",
                           highlightthickness=0)
        canvas.pack(pady=4)

        decision = {cid: None for cid in new_masks}
        new_mips_small = {}
        for cid, m in new_masks.items():
            nmm = np.max(m.astype(np.uint8), axis=0)
            small = np.array(
                Image.fromarray(nmm * 255, mode="L")
                .resize((disp_w, disp_h), Image.NEAREST)) > 127
            new_mips_small[cid] = small

        colors = {}
        for i, cid in enumerate(new_masks.keys()):
            hue = (i * 0.618033988749895) % 1.0
            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
            colors[cid] = (r * 255, g * 255, b * 255)

        photo_ref = {"img": None}

        def render():
            rgb = np.stack([base_small, base_small, base_small], axis=-1).astype(np.float32)
            for cid, d in decision.items():
                mask = new_mips_small[cid]
                if not mask.any():
                    continue
                if d == "keep":
                    col = np.array([40, 255, 40], dtype=np.float32); a = 0.65
                elif d == "delete":
                    col = np.array([255, 40, 40], dtype=np.float32); a = 0.7
                else:
                    col = np.array(colors[cid], dtype=np.float32); a = 0.45
                rgb[mask] = rgb[mask] * (1 - a) + col * a
            pil = Image.fromarray(rgb.astype(np.uint8), mode="RGB")
            photo = ImageTk.PhotoImage(pil)
            photo_ref["img"] = photo
            canvas.delete("all")
            canvas.create_image(0, 0, anchor="nw", image=photo)

        def on_click(event):
            x = int(event.x / scale)
            y = int(event.y / scale)
            if not (0 <= x < x_dim and 0 <= y < y_dim):
                return
            cid = None
            for c, m in new_masks.items():
                if m[:, y, x].any():
                    cid = c
                    break
            if cid is None:
                return
            if event.num == 1:
                decision[cid] = None if decision[cid] == "keep" else "keep"
            elif event.num == 3:
                decision[cid] = None if decision[cid] == "delete" else "delete"
            render()

        canvas.bind("<Button-1>", on_click)
        canvas.bind("<Button-3>", on_click)
        render()

        def finish():
            self.app.stop_work_panel()
            self._ok(decision)

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            help_text=self.app._t("dlg.review.help")
        )

    def _dialog_soma_repick(self, kw):
        mip = kw["mip"]; old_mip = kw.get("old_mip")
        panel = self.app.start_work_panel(self.app._t("dlg.repick.title"))

        arr = self._to_uint8(mip)
        h, w = arr.shape
        disp_w, disp_h, scale = self._fit_size(w, h, max_w=500, max_h=470)

        base_pil = self._pil_from_gray(arr).resize((disp_w, disp_h), Image.BILINEAR)
        base_photo = ImageTk.PhotoImage(base_pil)

        canvas = tk.Canvas(panel, width=disp_w, height=disp_h, bg="black",
                           highlightthickness=0)
        canvas.pack(pady=4)
        canvas.create_image(0, 0, anchor="nw", image=base_photo)
        canvas.image = base_photo

        pt = {"val": None}
        marker = {"id": None}

        def draw_old_outline():
            if old_mip is None:
                return
            small = np.array(
                Image.fromarray((old_mip > 0).astype(np.uint8) * 255, mode="L")
                .resize((disp_w, disp_h), Image.NEAREST)) > 127
            ys, xs = np.where(small)
            if len(xs) == 0:
                return
            step = max(1, len(xs) // 400)
            for i in range(0, len(xs), step):
                canvas.create_rectangle(xs[i], ys[i], xs[i] + 1, ys[i] + 1,
                                        outline="red")

        draw_old_outline()

        def on_click(event):
            x = int(event.x / scale)
            y = int(event.y / scale)
            pt["val"] = (y, x)
            if marker["id"] is not None:
                canvas.delete(marker["id"])
            r = 6
            marker["id"] = canvas.create_oval(event.x - r, event.y - r,
                                              event.x + r, event.y + r,
                                              outline="white", fill="lime", width=2)

        canvas.bind("<Button-1>", on_click)

        def finish():
            if pt["val"] is None:
                messagebox.showwarning(self.app._t("dlg.repick.title"),
                                       self.app._t("dlg.repick.point_warning"))
                return
            self.app.stop_work_panel()
            self._ok(pt["val"])

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.done"),
            cancel=True, cancel_value=None,
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.repick.help")
        )

class RoundedButton(tk.Canvas):
    def __init__(self, parent, text="", command=None,
                 width=160, height=36, radius=18,
                 bg="#7aa2f7", fg="#1a1c22",
                 hover_bg="#5d84d6", hover_fg=None,
                 active_bg="#4a6fbf",
                 font=("Segoe UI", 10, "bold"),
                 pad_x=16, disabled_bg="#3d4048",
                 disabled_fg="#6a6e78"):
        parent_bg = parent.cget("bg") if "bg" in parent.keys() else "#23252b"
        super().__init__(parent, width=width, height=height,
                         bg=parent_bg, highlightthickness=0, bd=0)
        self._label = text
        self._command = command
        self._btn_w = width
        self._btn_h = height
        self._btn_r = radius
        self._bg = bg
        self._fg = fg
        self._hover_bg = hover_bg
        self._hover_fg = hover_fg if hover_fg is not None else fg
        self._active_bg = active_bg
        self._disabled_bg = disabled_bg
        self._disabled_fg = disabled_fg
        self._font = font
        self._pad_x = pad_x
        self._state = "normal"
        self._is_hover = False
        self._is_press = False

        self._draw()
        self.bind("<Enter>", self._on_enter)
        self.bind("<Leave>", self._on_leave)
        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<ButtonRelease-1>", self._on_release)

    def _rounded_rect(self, x1, y1, x2, y2, r, **kwargs):
        pts = [
            x1 + r, y1,
            x2 - r, y1,
            x2, y1,
            x2, y1 + r,
            x2, y2 - r,
            x2, y2,
            x2 - r, y2,
            x1 + r, y2,
            x1, y2,
            x1, y2 - r,
            x1, y1 + r,
            x1, y1,
            ]
        return self.create_polygon(pts, smooth=True, **kwargs)

    def _draw(self):
        self.delete("all")
        if self._state == "disabled":
            bg = self._disabled_bg
            fg = self._disabled_fg
        elif self._is_press:
            bg = self._active_bg
            fg = self._hover_fg
        elif self._is_hover:
            bg = self._hover_bg
            fg = self._hover_fg
        else:
            bg = self._bg
            fg = self._fg

        self._rounded_rect(1, 1, self._btn_w - 1, self._btn_h - 1,
                           self._btn_r, fill=bg, outline="")
        self.create_text(self._btn_w // 2, self._btn_h // 2,
                         text=self._label, fill=fg, font=self._font)

    def _on_enter(self, _):
        if self._state == "disabled":
            return
        self._is_hover = True
        self.configure(cursor="hand2")
        self._draw()

    def _on_leave(self, _):
        self._is_hover = False
        self._is_press = False
        self.configure(cursor="")
        self._draw()

    def _on_press(self, _):
        if self._state == "disabled":
            return
        self._is_press = True
        self._draw()

    def _on_release(self, event):
        if self._state == "disabled":
            return
        was_press = self._is_press
        self._is_press = False
        self._draw()
        if was_press and 0 <= event.x <= self._btn_w and 0 <= event.y <= self._btn_h:
            if self._command:
                self._command()

    def config_state(self, state):
        self._state = state
        self._draw()

    def set_text(self, text):
        self._label = text
        self._draw()


class RoundedFrame(tk.Frame):
    def __init__(self, parent, radius=14, bg="#2b2e36",
                 border="#3d4048", inner_pad=10, **kwargs):
        parent_bg = parent.cget("bg") if "bg" in parent.keys() else "#23252b"
        super().__init__(parent, bg=parent_bg, bd=0,
                         highlightthickness=0, **kwargs)
        self._radius = radius
        self._card_bg = bg
        self._border = border
        self._inner_pad = inner_pad

        self._canvas = tk.Canvas(self, bg=parent_bg,
                                 highlightthickness=0, bd=0)
        self._canvas.place(x=0, y=0, relwidth=1, relheight=1)

        self._inner = tk.Frame(self, bg=bg)
        self._inner.pack(fill="both", expand=True,
                         padx=inner_pad, pady=inner_pad)

        self.bind("<Configure>", lambda e: self._redraw())

    def _redraw(self):
        self._canvas.delete("all")
        w = self.winfo_width()
        h = self.winfo_height()
        if w <= 2 or h <= 2:
            return
        r = self._radius
        pts = [
            1 + r, 1, w - 1 - r, 1, w - 1, 1,
            w - 1, 1 + r, w - 1, h - 1 - r, w - 1, h - 1,
            w - 1 - r, h - 1, 1 + r, h - 1, 1, h - 1,
            1, h - 1 - r, 1, 1 + r, 1, 1,
            ]
        self._canvas.create_polygon(
            pts, smooth=True,
            fill=self._card_bg, outline=self._border, width=1,
        )
        self._inner.configure(bg=self._card_bg)

    @property
    def inner(self):
        return self._inner

    def bind_configure(self):
        self.bind("<Configure>", lambda e: self._redraw())

class Chip(tk.Canvas):
    def __init__(self, parent, text="", value=None, variable=None,
                 bg="#33363f", bg_active="#7aa2f7",
                 fg="#d4d4dc", fg_active="#1a1c22",
                 font=("Segoe UI", 9), padx=12, pady=6):
        parent_bg = parent.cget("bg") if "bg" in parent.keys() else "#23252b"
        super().__init__(parent, bg=parent_bg, highlightthickness=0, bd=0)
        self._label = text
        self._value = value
        self._var = variable
        self._bg = bg
        self._bg_active = bg_active
        self._fg = fg
        self._fg_active = fg_active
        self._font = font
        self._padx = padx
        self._pady = pady
        self._hover = False

        self._measure()
        self._draw()
        self.bind("<Enter>", self._on_enter)
        self.bind("<Leave>", self._on_leave)
        self.bind("<Button-1>", self._on_click)
        if self._var is not None:
            self._var.trace_add("write", lambda *_: self._draw())

    def _measure(self):
        tmp = self.create_text(0, 0, text=self._label, font=self._font, anchor="nw")
        bbox = self.bbox(tmp)
        self.delete(tmp)
        w = (bbox[2] - bbox[0]) + self._padx * 2
        h = (bbox[3] - bbox[1]) + self._pady * 2
        self.configure(width=w, height=h)
        self._chip_w = w
        self._chip_h = h

    def _rounded_rect(self, x1, y1, x2, y2, r, **kw):
        pts = [
            x1 + r, y1, x2 - r, y1, x2, y1,
            x2, y1 + r, x2, y2 - r, x2, y2,
            x2 - r, y2, x1 + r, y2, x1, y2,
            x1, y2 - r, x1, y1 + r, x1, y1,
            ]
        return self.create_polygon(pts, smooth=True, **kw)

    def _is_selected(self):
        if self._var is None:
            return False
        return str(self._var.get()) == str(self._value)

    def _draw(self):
        self.delete("all")
        selected = self._is_selected()
        if selected:
            bg = self._bg_active
            fg = self._fg_active
        elif self._hover:
            bg = "#3c404b"
            fg = self._fg
        else:
            bg = self._bg
            fg = self._fg
        r = self._chip_h // 2
        self._rounded_rect(1, 1, self._chip_w - 1, self._chip_h - 1, r,
                           fill=bg, outline="")
        self.create_text(self._chip_w // 2, self._chip_h // 2,
                         text=self._label, fill=fg, font=self._font)

    def _on_enter(self, _):
        self._hover = True
        self.configure(cursor="hand2")
        self._draw()

    def _on_leave(self, _):
        self._hover = False
        self.configure(cursor="")
        self._draw()

    def _on_click(self, _):
        if self._var is not None:
            self._var.set(self._value)

    def set_label(self, text):
        self._label = text
        self._measure()
        self._draw()

class RoundedScrollbar(tk.Canvas):
    def __init__(self, parent, command=None, width=10, orient="vertical",
                 bg="#23252b", thumb="#4a4d57", thumb_hover="#7aa2f7", **kwargs):
        if orient == "horizontal":
            kwargs["height"] = width
            super().__init__(parent, bg=bg, highlightthickness=0, bd=0, **kwargs)
        else:
            kwargs["width"] = width
            super().__init__(parent, bg=bg, highlightthickness=0, bd=0, **kwargs)

        self._command = command
        self._orient = orient
        self._bg = bg
        self._thumb = thumb
        self._thumb_hover = thumb_hover
        self._first = 0.0
        self._last = 1.0
        self._hover = False
        self._drag = False

        self.bind("<Configure>", lambda e: self._redraw())
        self.bind("<Enter>", self._on_enter)
        self.bind("<Leave>", self._on_leave)
        self.bind("<ButtonPress-1>", self._on_press)
        self.bind("<B1-Motion>", self._on_drag)
        self.bind("<ButtonRelease-1>", self._on_release)

    def set(self, first, last):
        try:
            self._first = float(first)
            self._last = float(last)
        except Exception:
            return
        self._redraw()

    def _rounded(self, x1, y1, x2, y2, r, **kw):
        pts = [
            x1 + r, y1, x2 - r, y1, x2, y1,
            x2, y1 + r, x2, y2 - r, x2, y2,
            x2 - r, y2, x1 + r, y2, x1, y2,
            x1, y2 - r, x1, y1 + r, x1, y1,
            ]
        return self.create_polygon(pts, smooth=True, **kw)

    def _redraw(self):
        self.delete("all")
        w = self.winfo_width()
        h = self.winfo_height()
        if w <= 1 or h <= 1:
            return

        if self._orient == "horizontal":
            left = int(self._first * w)
            right = int(self._last * w)
            if right - left < 24:
                mid = (left + right) // 2
                left = max(0, mid - 12)
                right = min(w, mid + 12)
            color = self._thumb_hover if (self._hover or self._drag) else self._thumb
            self._rounded(left + 1, 1, right - 1, h - 1, (h - 2) // 2,
                          fill=color, outline="")
        else:
            top = int(self._first * h)
            bot = int(self._last * h)
            if bot - top < 24:
                mid = (top + bot) // 2
                top = max(0, mid - 12)
                bot = min(h, mid + 12)
            color = self._thumb_hover if (self._hover or self._drag) else self._thumb
            self._rounded(1, top + 1, w - 1, bot - 1, (w - 2) // 2,
                          fill=color, outline="")

    def _on_enter(self, _):
        self._hover = True
        self._redraw()

    def _on_leave(self, _):
        self._hover = False
        self._redraw()

    def _on_press(self, event):
        self._drag = True
        self._scroll_to(event.x if self._orient == "horizontal" else event.y)
        self._redraw()

    def _on_drag(self, event):
        if self._drag:
            self._scroll_to(event.x if self._orient == "horizontal" else event.y)

    def _on_release(self, _):
        self._drag = False
        self._redraw()

    def _scroll_to(self, val):
        if self._command is None:
            return
        length = self.winfo_width() if self._orient == "horizontal" else self.winfo_height()
        if length <= 1:
            return
        frac = max(0.0, min(1.0, val / length))
        self._command("moveto", str(frac))

class RoundedSlider(tk.Canvas):
    def __init__(self, parent, from_=0, to=100, variable=None,
                 length=400, height=36, command=None,
                 bg="#33363f", track="#3d4048", fill="#7aa2f7",
                 thumb="#d4d4dc", thumb_hover="#ffffff",
                 showvalue=True, font=("Segoe UI", 9, "bold"),
                 fg="#d4d4dc"):
        super().__init__(parent, width=length, height=height,
                         bg=bg, highlightthickness=0, bd=0)
        self._from = from_
        self._to = to
        self._var = variable
        self._cmd = command
        self._len = length
        self._h = height
        self._bg = bg
        self._track = track
        self._fill = fill
        self._thumb = thumb
        self._thumb_hover = thumb_hover
        self._showvalue = showvalue
        self._font = font
        self._fg = fg
        self._hover = False
        self._drag = False

        if self._var is None:
            self._var = tk.DoubleVar(value=from_)
        self._var.trace_add("write", lambda *_: self._redraw())

        self.bind("<Configure>", lambda e: self._redraw())
        self.bind("<Button-1>", self._on_press)
        self.bind("<B1-Motion>", self._on_drag)
        self.bind("<ButtonRelease-1>", self._on_release)
        self.bind("<Enter>", lambda e: (setattr(self, "_hover", True), self._redraw()))
        self.bind("<Leave>", lambda e: (setattr(self, "_hover", False), self._redraw()))

    def _rounded(self, x1, y1, x2, y2, r, **kw):
        pts = [
            x1 + r, y1, x2 - r, y1, x2, y1,
            x2, y1 + r, x2, y2 - r, x2, y2,
            x2 - r, y2, x1 + r, y2, x1, y2,
            x1, y2 - r, x1, y1 + r, x1, y1,
            ]
        return self.create_polygon(pts, smooth=True, **kw)

    def _frac(self):
        try:
            v = float(self._var.get())
        except Exception:
            v = self._from
        if self._to == self._from:
            return 0.0
        return max(0.0, min(1.0, (v - self._from) / (self._to - self._from)))

    def _set_from_x(self, x):
        pad = 14
        w = self._len
        frac = (x - pad) / max(1, (w - 2 * pad))
        frac = max(0.0, min(1.0, frac))
        v = self._from + frac * (self._to - self._from)
        self._var.set(int(round(v)) if self._to - self._from > 5 else v)
        if self._cmd:
            self._cmd(self._var.get())

    def _redraw(self):
        self.delete("all")
        w = self._len
        h = self._h
        cy = h // 2
        pad = 14
        track_h = 6
        r = track_h // 2

        self._rounded(pad, cy - r, w - pad, cy + r, r,
                      fill=self._track, outline="")

        f = self._frac()
        x_thumb = pad + f * (w - 2 * pad)
        if x_thumb > pad + 2:
            self._rounded(pad, cy - r, int(x_thumb), cy + r, r,
                          fill=self._fill, outline="")

        thumb_r = 10
        color = self._thumb_hover if (self._hover or self._drag) else self._thumb
        self.create_oval(x_thumb - thumb_r, cy - thumb_r,
                         x_thumb + thumb_r, cy + thumb_r,
                         fill=color, outline=self._bg, width=2)

        if self._showvalue:
            self.create_text(w - pad, cy - thumb_r - 6, anchor="e",
                             text=str(self._var.get()),
                             fill=self._fg, font=self._font)

    def _on_press(self, event):
        self._drag = True
        self._set_from_x(event.x)
        self._redraw()

    def _on_drag(self, event):
        if self._drag:
            self._set_from_x(event.x)
            self._redraw()

    def _on_release(self, _):
        self._drag = False
        self._redraw()


class StyledDialog(tk.Toplevel):
    def __init__(self, parent, title, text, kind="info",
                 buttons=(("OK", "ok"),), bg="#23252b",
                 fg="#d4d4dc", fg_dim="#8e93a0",
                 accent="#7aa2f7", err="#e07a7a",
                 warn="#e6b673", border="#3d4048"):
        super().__init__(parent)
        self.title(title)
        self.configure(bg=bg)
        self.resizable(False, False)
        self.transient(parent)
        self.grab_set()
        self.result = "ok"

        self._bg = bg
        self._fg = fg
        self._fg_dim = fg_dim
        self._accent = accent
        self._err = err
        self._warn = warn
        self._border = border
        self._buttons = buttons

        outer = tk.Frame(self, bg=bg, bd=0,
                         highlightthickness=1,
                         highlightbackground=border,
                         highlightcolor=border)
        outer.pack(fill="both", expand=True, padx=10, pady=10)

        icon_color = {
            "info": accent,
            "warn": warn,
            "error": err,
        }.get(kind, accent)

        header = tk.Canvas(outer, width=40, height=40, bg=bg,
                           highlightthickness=0, bd=0)
        header.pack(side="left", padx=(12, 8), pady=14, anchor="n")
        header.create_oval(4, 4, 36, 36, outline=icon_color, width=2)
        header.create_text(20, 20,
                           text="!" if kind == "warn" else ("×" if kind == "error" else "i"),
                           fill=icon_color, font=("Segoe UI", 14, "bold"))

        body = tk.Frame(outer, bg=bg)
        body.pack(side="left", fill="both", expand=True,
                  padx=(0, 14), pady=(14, 4))

        tk.Label(body, text=title, bg=bg, fg=fg,
                 font=("Segoe UI", 11, "bold"), anchor="w").pack(fill="x")
        tk.Label(body, text=text, bg=bg, fg=fg_dim,
                 font=("Segoe UI", 10), justify="left",
                 wraplength=360, anchor="w").pack(fill="x", pady=(6, 0))

        btns = tk.Frame(body, bg=bg)
        btns.pack(fill="x", pady=(14, 0))

        def _mk(text, val, primary=False, danger=False):
            if danger:
                b, f, hb = "#3b2a2e", err, "#4d3238"
            elif primary:
                b, f, hb = accent, "#1a1c22", "#5d84d6"
            else:
                b, f, hb = "#33363f", fg, "#3c404b"
            btn = RoundedButton(btns, text=text,
                                command=lambda: self._choose(val),
                                width=120, height=32, radius=16,
                                bg=b, fg=f,
                                hover_bg=hb, hover_fg=f,
                                active_bg=hb,
                                font=("Segoe UI", 9,
                                      "bold" if primary or danger else "normal"))
            btn.pack(side="right", padx=4)
            return btn

        for i, (text, val) in enumerate(buttons):
            _mk(text, val, primary=(i == 0 and kind != "error"),
                danger=(kind == "error" and i == 0))

        self.protocol("WM_DELETE_WINDOW", lambda: self._choose("cancel"))
        self.update_idletasks()
        pw = parent.winfo_rootx() + parent.winfo_width() // 2
        ph = parent.winfo_rooty() + parent.winfo_height() // 2
        w = self.winfo_width()
        h = self.winfo_height()
        self.geometry(f"+{pw - w // 2}+{ph - h // 2}")

    def _choose(self, value):
        self.result = value
        self.destroy()

    @classmethod
    def show(cls, parent, title, text, kind="info", buttons=(("OK", "ok"),)):
        dlg = cls(parent, title, text, kind=kind, buttons=buttons)
        parent.wait_window(dlg)
        return dlg.result

class App:
    BG        = "#23252b"
    BG_SOFT   = "#2b2e36"
    BG_CARD   = "#33363f"
    BG_HOVER  = "#3c404b"

    FG        = "#d4d4dc"
    FG_DIM    = "#8e93a0"
    FG_MUTED  = "#6a6e78"

    ACCENT    = "#7aa2f7"
    ACCENT_D  = "#5d84d6"
    OK        = "#8fd694"
    WARN      = "#e6b673"
    ERR       = "#e07a7a"
    BORDER    = "#3d4048"

    def __init__(self):
        self.root = tk.Tk()
        self.root.title("ASTRO-ANALYZER")
        self.root.minsize(1180, 720)
        self.root.configure(bg=self.BG)

        try:
            self.root.state("zoomed")
        except Exception:
            sw = self.root.winfo_screenwidth()
            sh = self.root.winfo_screenheight()
            self.root.geometry(f"{sw}x{sh}+0+0")

        self._init_style()

        self.bridge = GuiBridge(self)
        self.worker = None

        self._current_lang = "ru"

        self._settings_labels = []
        self._settings_buttons = []
        self._settings_checks = []
        self._settings_chips = []
        self._settings_static_labels = []
        self._src_labels = []
        self._log_title_label = None

        body = tk.Frame(self.root, bg=self.BG)
        body.pack(fill="both", expand=True, padx=18, pady=(14, 14))

        self._build_workspace(body)

        self._work_panel_active = False
        self._analysis_mode = False

    def _t(self, key, **kwargs):
        lang = getattr(self, "var_lang", None)
        cur = lang.get() if lang is not None else "ru"
        table = TRANSLATIONS.get(cur, TRANSLATIONS["ru"])
        template = table.get(key)
        if template is None:
            template = TRANSLATIONS["ru"].get(key, key)
        try:
            return template.format(**kwargs)
        except Exception:
            return template

    def _apply_language(self, lang):
        self._current_lang = lang
        self.root.title(self._t("app.title"))

        if hasattr(self, "_header_canvas"):
            try:
                self._header_canvas.event_generate("<Configure>")
            except Exception:
                pass

        for lbl, key in getattr(self, "_settings_labels", []):
            try:
                lbl.config(text=self._t(key))
            except Exception:
                pass

        for btn, key in getattr(self, "_settings_buttons", []):
            try:
                btn.set_text(self._t(key))
            except Exception:
                pass

        for chk, key in getattr(self, "_settings_checks", []):
            try:
                chk.config(text=self._t(key))
            except Exception:
                pass

        for chip, key in getattr(self, "_settings_chips", []):
            try:
                chip.set_label(self._t(key))
            except Exception:
                pass

        for lbl, key in getattr(self, "_settings_static_labels", []):
            try:
                lbl.config(text=self._t(key))
            except Exception:
                pass

        for lbl, key in getattr(self, "_src_labels", []):
            try:
                cur = lbl.cget("text")
                prev_ru = TRANSLATIONS["ru"].get(key, cur)
                prev_en = TRANSLATIONS["en"].get(key, cur)
                if cur in (prev_ru, prev_en):
                    lbl.config(text=self._t(key))
            except Exception:
                pass

        if hasattr(self, "btn_run"):
            try:
                self.btn_run.set_text(self._t("ui.settings.run"))
            except Exception:
                pass

        if hasattr(self, "btn_exit"):
            try:
                self.btn_exit.set_text(self._t("ui.settings.exit"))
            except Exception:
                pass

        if hasattr(self, "btn_back_welcome"):
            try:
                self.btn_back_welcome.set_text(self._t("ui.settings.back_to_settings"))
            except Exception:
                pass

        if hasattr(self, "btn_to_log"):
            try:
                self.btn_to_log.set_text(self._t("ui.settings.back_to_log"))
            except Exception:
                pass

        if getattr(self, "_log_title_label", None) is not None:
            try:
                self._log_title_label.config(text=self._t("ui.settings.log_title"))
            except Exception:
                pass

    def _init_style(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure("TFrame", background=self.BG)
        style.configure("TLabel", background=self.BG, foreground=self.FG,
                        font=("Segoe UI", 10))
        style.configure("TRadiobutton", background=self.BG, foreground=self.FG,
                        font=("Segoe UI", 10))
        style.map("TRadiobutton",
                  background=[("active", self.BG)],
                  foreground=[("active", self.FG)])
        style.configure("TCheckbutton", background=self.BG, foreground=self.FG,
                        font=("Segoe UI", 10))
        style.map("TCheckbutton",
                  background=[("active", self.BG)],
                  foreground=[("active", self.FG)])
        style.configure("TEntry", fieldbackground=self.BG_SOFT,
                        foreground=self.FG, insertcolor=self.FG,
                        bordercolor=self.BORDER, lightcolor=self.BORDER,
                        darkcolor=self.BORDER, borderwidth=1)
        style.configure("TScale", background=self.BG_CARD,
                        troughcolor=self.BG_SOFT, bordercolor=self.BG_SOFT,
                        lightcolor=self.ACCENT, darkcolor=self.ACCENT)
        style.configure("TLabelframe", background=self.BG, bordercolor=self.BORDER,
                        relief="solid", borderwidth=1)
        style.configure("TLabelframe.Label", background=self.BG, foreground=self.FG,
                        font=("Segoe UI", 10, "bold"))

    def _rounded_rect_canvas(self, canvas, x1, y1, x2, y2, r, **kw):
        pts = [
            x1 + r, y1, x2 - r, y1, x2, y1,
            x2, y1 + r, x2, y2 - r, x2, y2,
            x2 - r, y2, x1 + r, y2, x1, y2,
            x1, y2 - r, x1, y1 + r, x1, y1,
            ]
        return canvas.create_polygon(pts, smooth=True, **kw)

    def _build_workspace(self, parent):
        self.right = tk.Frame(parent, width=600, bg=self.BG_CARD,
                              bd=0, highlightthickness=1,
                              highlightbackground=self.BORDER,
                              highlightcolor=self.BORDER)
        self.right.pack_propagate(False)
        self.right.pack(side="right", fill="y")

        self.left_col = tk.Frame(parent, bg=self.BG)
        self.left_col.pack(side="left", fill="both", expand=True,
                           padx=(0, 14))

        self._build_header(self.left_col)

        self.left_content = tk.Frame(self.left_col, bg=self.BG)
        self.left_content.pack(fill="both", expand=True)

        self._build_welcome(self.left_content)
        self._build_log(self.left_content)
        self.log_container.pack_forget()

        self._build_footer(self.left_col)

    def _build_header(self, parent):
        header_row = tk.Frame(parent, bg=self.BG)
        header_row.pack(fill="x", pady=(0, 10))

        lang_frame = tk.Frame(header_row, bg=self.BG)
        lang_frame.pack(side="right", padx=(8, 0))

        lbl_lang = tk.Label(lang_frame, text=self._t("ui.header.lang"),
                            bg=self.BG, fg=self.FG_DIM,
                            font=("Segoe UI", 9))
        lbl_lang.pack(side="left", padx=(0, 6))
        self._settings_static_labels.append((lbl_lang, "ui.header.lang"))

        self.var_lang = tk.StringVar(value="ru")
        Chip(lang_frame, "RU", "ru", self.var_lang,
             bg=self.BG_CARD, bg_active=self.ACCENT).pack(side="left")
        Chip(lang_frame, "EN", "en", self.var_lang,
             bg=self.BG_CARD, bg_active=self.ACCENT).pack(side="left", padx=(6, 0))

        def _on_lang_change(*_):
            self._apply_language(self.var_lang.get())

        self.var_lang.trace_add("write", _on_lang_change)

        header = tk.Canvas(header_row, height=64, bg=self.BG,
                           highlightthickness=0, bd=0)
        header.pack(side="left", fill="x", expand=True)

        self._logo_photo = None
        logo_path_candidates = [
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "lab_logo.png"),
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "logo.png"),
        ]
        logo_path = next((p for p in logo_path_candidates if os.path.exists(p)), None)

        if logo_path is not None:
            try:
                pil = Image.open(logo_path).convert("RGBA")
                max_h = 44
                if pil.height > max_h:
                    ratio = max_h / pil.height
                    pil = pil.resize((int(pil.width * ratio), max_h),
                                     Image.LANCZOS)
                self._logo_photo = ImageTk.PhotoImage(pil)
            except Exception:
                self._logo_photo = None

        def draw(event=None):
            header.delete("all")
            w = header.winfo_width()
            if w <= 4:
                return
            self._rounded_rect_canvas(
                header, 0, 2, w, 62, 20,
                fill=self.BG_CARD, outline=self.BORDER, width=1
            )
            header.create_text(
                24, 32, anchor="w",
                text=self._t("app.header_title"),
                fill=self.FG, font=("Segoe UI", 14, "bold")
            )
            if self._logo_photo is not None:
                header.create_image(w - 24, 32, anchor="e",
                                    image=self._logo_photo)
            else:
                header.create_text(
                    w - 24, 32, anchor="e",
                    text="LAB",
                    fill=self.FG_DIM,
                    font=("Segoe UI", 13, "italic", "bold")
                )

        self._header_canvas = header
        header.bind("<Configure>", draw)

    def _build_footer(self, parent):
        footer = tk.Frame(parent, bg=self.BG)
        footer.pack(side="bottom", fill="x", pady=(10, 0))

        self.btn_run = RoundedButton(
            footer, self._t("ui.settings.run"),
            command=self.start,
            width=220, height=44, radius=22,
            bg=self.ACCENT, fg="#1a1c22",
            hover_bg=self.ACCENT_D, hover_fg="#ffffff",
            active_bg="#4a6fbf",
            font=("Segoe UI", 10, "bold")
        )
        self.btn_run.pack(side="left")

        self.btn_exit = RoundedButton(
            footer, self._t("ui.settings.exit"),
            command=self.root.quit,
            width=120, height=44, radius=22,
            bg=self.BG_CARD, fg=self.FG,
            hover_bg=self.BG_HOVER, hover_fg=self.FG,
            active_bg=self.BORDER,
            font=("Segoe UI", 10)
        )
        self.btn_exit.pack(side="right")

        self._settings_buttons.append((self.btn_run, "ui.settings.run"))
        self._settings_buttons.append((self.btn_exit, "ui.settings.exit"))

    def _build_welcome(self, parent):
        self.welcome_container = tk.Frame(parent, bg=self.BG)
        self.welcome_container.pack(fill="both", expand=True)

        wrap = tk.Frame(self.welcome_container, bg=self.BG_CARD,
                        bd=0, highlightthickness=1,
                        highlightbackground=self.BORDER,
                        highlightcolor=self.BORDER)
        wrap.pack(fill="both", expand=True)

        header = tk.Frame(wrap, bg=self.BG_CARD)
        header.pack(fill="x", padx=16, pady=(16, 4))

        self.btn_to_log = RoundedButton(
            header, self._t("ui.settings.back_to_log"),
            command=self._switch_to_analysis_mode,
            width=180, height=26, radius=13,
            bg=self.BG_CARD, fg=self.FG_DIM,
            hover_bg=self.BG_HOVER, hover_fg=self.FG,
            active_bg=self.BORDER,
            font=("Segoe UI", 9),
        )
        self.btn_to_log.pack_forget()

        lbl_main = tk.Label(header,
                            text=self._t("ui.settings.title_main"),
                            bg=self.BG_CARD, fg=self.FG,
                            font=("Segoe UI", 10, "bold"),
                            anchor="w")
        lbl_main.pack(fill="x")

        lbl_sub = tk.Label(header,
                           text=self._t("ui.settings.subtitle"),
                           bg=self.BG_CARD, fg=self.FG_DIM,
                           font=("Segoe UI", 9),
                           anchor="w", justify="left",
                           wraplength=900)
        lbl_sub.pack(fill="x", pady=(4, 0))

        self._settings_static_labels.append((lbl_main, "ui.settings.title_main"))
        self._settings_static_labels.append((lbl_sub, "ui.settings.subtitle"))

        body = tk.Frame(wrap, bg=self.BG_CARD)
        body.pack(fill="x", padx=16, pady=(8, 16))

        self._build_settings(body)

    def _build_log(self, parent):
        self.log_container = RoundedFrame(parent, radius=16, bg=self.BG_SOFT,
                                          border=self.BORDER, inner_pad=14)
        self.log_container.pack(fill="both", expand=True)

        title_row = tk.Frame(self.log_container.inner, bg=self.BG_SOFT)
        title_row.pack(fill="x", pady=(0, 8))

        self._log_title_label = tk.Label(title_row,
                                         text=self._t("ui.settings.log_title"),
                                         bg=self.BG_SOFT, fg=self.FG,
                                         font=("Segoe UI", 10, "bold"),
                                         anchor="w")
        self._log_title_label.pack(side="left")

        self.btn_back_welcome = RoundedButton(
            title_row, self._t("ui.settings.back_to_settings"),
            command=self._switch_to_welcome_mode,
            width=140, height=26, radius=13,
            bg=self.BG_CARD, fg=self.FG_DIM,
            hover_bg=self.BG_HOVER, hover_fg=self.FG,
            active_bg=self.BORDER,
            font=("Segoe UI", 9),
        )

        text_wrap = tk.Frame(self.log_container.inner, bg=self.BG_SOFT)
        text_wrap.pack(fill="both", expand=True)

        self.log_text = tk.Text(text_wrap, wrap="word",
                                state="normal",
                                bg=self.BG_SOFT, fg=self.FG,
                                insertbackground=self.ACCENT,
                                selectbackground=self.ACCENT_D,
                                selectforeground="#ffffff",
                                relief="flat", bd=0,
                                padx=4, pady=2,
                                font=("Consolas", 10))
        log_scroll = RoundedScrollbar(text_wrap, command=self.log_text.yview,
                                      width=10,
                                      bg=self.BG_SOFT,
                                      thumb=self.BORDER,
                                      thumb_hover=self.ACCENT)
        self.log_text.configure(yscrollcommand=log_scroll.set)
        self.log_text.pack(side="left", fill="both", expand=True)
        log_scroll.pack(side="right", fill="y")

        self.log_text.tag_configure("ok", foreground=self.OK)
        self.log_text.tag_configure("warn", foreground=self.WARN)
        self.log_text.tag_configure("err", foreground=self.ERR)
        self.log_text.tag_configure("accent", foreground=self.ACCENT)
        self.log_text.tag_configure("dim", foreground=self.FG_DIM)

    def _switch_to_analysis_mode(self):
        if self._analysis_mode:
            return
        self._analysis_mode = True
        self.welcome_container.pack_forget()
        self.log_container.pack(fill="both", expand=True)
        self.btn_back_welcome.pack(side="right")
        self.btn_to_log.pack_forget()

    def _switch_to_welcome_mode(self):
        if not self._analysis_mode:
            return
        self._analysis_mode = False
        self.btn_back_welcome.pack_forget()
        self.log_container.pack_forget()
        self.welcome_container.pack(fill="both", expand=True)
        self.btn_to_log.pack(side="right")

    def _build_settings(self, parent):
        s = tk.Frame(parent, bg=self.BG_CARD)
        s.pack(fill="x")

        row = self._row(s, "ui.settings.file")
        self.btn_browse = RoundedButton(row, self._t("ui.settings.browse"),
                                        command=self._choose_file,
                                        width=110, height=30, radius=15,
                                        bg=self.BG_HOVER, fg=self.FG,
                                        hover_bg=self.ACCENT, hover_fg="#1a1c22",
                                        active_bg=self.ACCENT_D,
                                        font=("Segoe UI", 9))
        self.btn_browse.pack(side="left", padx=(0, 8))
        self._settings_buttons.append((self.btn_browse, "ui.settings.browse"))

        self.var_file_path = tk.StringVar(value="")
        self._entry(row, self.var_file_path, width=10).pack(
            side="left", fill="x", expand=True)

        row = self._row(s, "ui.settings.load_mode")
        self.var_load_mode = tk.StringVar(value="one")
        self.chip_load_one = Chip(row, self._t("ui.settings.load_one"), "one",
                                  self.var_load_mode)
        self.chip_load_one.pack(side="left")
        self.chip_load_dir = Chip(row, self._t("ui.settings.load_dir"), "dir",
                                  self.var_load_mode)
        self.chip_load_dir.pack(side="left", padx=8)
        self._settings_chips = getattr(self, "_settings_chips", [])
        self._settings_chips.append((self.chip_load_one, "ui.settings.load_one"))
        self._settings_chips.append((self.chip_load_dir, "ui.settings.load_dir"))

        row = self._row(s, "ui.settings.metadata")
        self.btn_read_metadata = RoundedButton(row, self._t("ui.settings.read_metadata"),
                                               command=self._read_metadata,
                                               width=190, height=30, radius=15,
                                               bg=self.BG_HOVER, fg=self.FG,
                                               hover_bg=self.ACCENT, hover_fg="#1a1c22",
                                               active_bg=self.ACCENT_D,
                                               font=("Segoe UI", 9))
        self.btn_read_metadata.pack(side="left")
        self._settings_buttons.append((self.btn_read_metadata, "ui.settings.read_metadata"))

        row = self._row(s, "ui.settings.cell_type")
        self.var_cell_type = tk.StringVar(value="astrocytes")
        self.chip_astro = Chip(row, self._t("cell.astrocytes_label"), "astrocytes",
                               self.var_cell_type)
        self.chip_astro.pack(side="left")
        self.chip_microglia = Chip(row, self._t("cell.microglia_label"), "microglia",
                                   self.var_cell_type)
        self.chip_microglia.pack(side="left", padx=8)
        self.chip_neurons = Chip(row, self._t("cell.neurons_label"), "neurons",
                                 self.var_cell_type)
        self.chip_neurons.pack(side="left")
        self._settings_chips.append((self.chip_astro, "cell.astrocytes_label"))
        self._settings_chips.append((self.chip_microglia, "cell.microglia_label"))
        self._settings_chips.append((self.chip_neurons, "cell.neurons_label"))

        row = self._row(s, "ui.settings.image_mode")
        self.var_image_mode = tk.StringVar(value="confocal")
        self.chip_confocal = Chip(row, self._t("mode.confocal"), "confocal",
                                  self.var_image_mode)
        self.chip_confocal.pack(side="left")
        self.chip_epi = Chip(row, self._t("mode.epifluorescence"), "epifluorescence",
                             self.var_image_mode)
        self.chip_epi.pack(side="left", padx=8)
        self._settings_chips.append((self.chip_confocal, "mode.confocal"))
        self._settings_chips.append((self.chip_epi, "mode.epifluorescence"))

        row = self._row(s, "ui.settings.deconv")
        self.var_use_deconv = tk.BooleanVar(value=True)
        chk = self._check(row, self._t("ui.settings.use_deconv"),
                          self.var_use_deconv)
        chk.pack(side="left")
        self._settings_checks.append((chk, "ui.settings.use_deconv"))

        row = self._row(s, "ui.settings.xy")
        self.var_xy = tk.StringVar(value="0.144")
        self.entry_xy = self._entry(row, self.var_xy, width=10)
        self.entry_xy.pack(side="left")
        self.lbl_xy_src = tk.Label(row, text=self._t("ui.settings.default"),
                                   bg=self.BG_CARD, fg=self.FG_DIM,
                                   font=("Segoe UI", 9, "italic"))
        self.lbl_xy_src.pack(side="left", padx=(8, 0))
        self._src_labels.append((self.lbl_xy_src, "ui.settings.default"))

        row = self._row(s, "ui.settings.z")
        self.var_z = tk.StringVar(value="1.0")
        self.entry_z = self._entry(row, self.var_z, width=10)
        self.entry_z.pack(side="left")
        self.lbl_z_src = tk.Label(row, text=self._t("ui.settings.default"),
                                  bg=self.BG_CARD, fg=self.FG_DIM,
                                  font=("Segoe UI", 9, "italic"))
        self.lbl_z_src.pack(side="left", padx=(8, 0))
        self._src_labels.append((self.lbl_z_src, "ui.settings.default"))

        hint = tk.Label(s, text=self._t("ui.settings.scale_hint"),
                        bg=self.BG_CARD, fg=self.FG_DIM,
                        font=("Segoe UI", 9, "italic"),
                        anchor="w", justify="left",
                        wraplength=760)
        hint.pack(fill="x", pady=(0, 6))
        self._settings_static_labels.append((hint, "ui.settings.scale_hint"))

        row = self._row(s, "ui.settings.cells", width=14)
        lbl_na_c = tk.Label(row, text=self._t("ui.settings.na"),
                            bg=self.BG_CARD, fg=self.FG,
                            font=("Segoe UI", 10))
        lbl_na_c.pack(side="left", padx=(0, 4))
        self._settings_labels.append((lbl_na_c, "ui.settings.na"))
        self.var_na_cells = tk.StringVar(value="0.95")
        self._entry(row, self.var_na_cells, width=6).pack(side="left")
        lbl_wl_c = tk.Label(row, text=self._t("ui.settings.lambda_nm"),
                            bg=self.BG_CARD, fg=self.FG,
                            font=("Segoe UI", 10))
        lbl_wl_c.pack(side="left", padx=(14, 4))
        self._settings_labels.append((lbl_wl_c, "ui.settings.lambda_nm"))
        self.var_wl_cells = tk.StringVar(value="699.5")
        self._entry(row, self.var_wl_cells, width=8).pack(side="left")
        lbl_eps_c = tk.Label(row, text=self._t("ui.settings.epsilon"),
                             bg=self.BG_CARD, fg=self.FG,
                             font=("Segoe UI", 10))
        lbl_eps_c.pack(side="left", padx=(14, 4))
        self._settings_labels.append((lbl_eps_c, "ui.settings.epsilon"))
        self.var_eps_cells = tk.StringVar(value="0.005")
        self._entry(row, self.var_eps_cells, width=8).pack(side="left")
        self.lbl_cells_src = tk.Label(row, text=self._t("ui.settings.default"),
                                      bg=self.BG_CARD, fg=self.FG_DIM,
                                      font=("Segoe UI", 9, "italic"))
        self.lbl_cells_src.pack(side="left", padx=(12, 0))
        self._src_labels.append((self.lbl_cells_src, "ui.settings.default"))

        row = self._row(s, "ui.settings.dapi", width=14)
        lbl_na_d = tk.Label(row, text=self._t("ui.settings.na"),
                            bg=self.BG_CARD, fg=self.FG,
                            font=("Segoe UI", 10))
        lbl_na_d.pack(side="left", padx=(0, 4))
        self._settings_labels.append((lbl_na_d, "ui.settings.na"))
        self.var_na_dapi = tk.StringVar(value="0.95")
        self._entry(row, self.var_na_dapi, width=6).pack(side="left")
        lbl_wl_d = tk.Label(row, text=self._t("ui.settings.lambda_nm"),
                            bg=self.BG_CARD, fg=self.FG,
                            font=("Segoe UI", 10))
        lbl_wl_d.pack(side="left", padx=(14, 4))
        self._settings_labels.append((lbl_wl_d, "ui.settings.lambda_nm"))
        self.var_wl_dapi = tk.StringVar(value="460")
        self._entry(row, self.var_wl_dapi, width=8).pack(side="left")
        lbl_eps_d = tk.Label(row, text=self._t("ui.settings.epsilon"),
                             bg=self.BG_CARD, fg=self.FG,
                             font=("Segoe UI", 10))
        lbl_eps_d.pack(side="left", padx=(14, 4))
        self._settings_labels.append((lbl_eps_d, "ui.settings.epsilon"))
        self.var_eps_dapi = tk.StringVar(value="0.005")
        self._entry(row, self.var_eps_dapi, width=8).pack(side="left")
        self.lbl_dapi_src = tk.Label(row, text=self._t("ui.settings.default"),
                                     bg=self.BG_CARD, fg=self.FG_DIM,
                                     font=("Segoe UI", 9, "italic"))
        self.lbl_dapi_src.pack(side="left", padx=(12, 0))
        self._src_labels.append((self.lbl_dapi_src, "ui.settings.default"))

        row = self._row(s, "ui.settings.filter")
        self.var_do_filter = tk.BooleanVar(value=True)
        chk_filter = self._check(row, self._t("ui.settings.run_filter"),
                                 self.var_do_filter)
        chk_filter.pack(side="left")
        self._settings_checks.append((chk_filter, "ui.settings.run_filter"))

        self._settings_widgets = []
        self._collect_settings_widgets(s)

    def _row(self, parent, label_key, width=26):
        row = tk.Frame(parent, bg=self.BG_CARD)
        row.pack(fill="x", pady=3)
        lbl = tk.Label(row, text=self._t(label_key), width=width, anchor="w",
                       bg=self.BG_CARD, fg=self.FG,
                       font=("Segoe UI", 10))
        lbl.pack(side="left")
        self._settings_labels.append((lbl, label_key))
        return row

    def _entry(self, parent, var, width=10):
        return tk.Entry(parent, textvariable=var, width=width,
                        bg=self.BG_SOFT, fg=self.FG,
                        insertbackground=self.ACCENT,
                        relief="flat", bd=0,
                        highlightthickness=1,
                        highlightbackground=self.BORDER,
                        highlightcolor=self.ACCENT,
                        font=("Consolas", 10))

    def _check(self, parent, text, var):
        return tk.Checkbutton(parent, text=text, variable=var,
                              bg=self.BG_CARD, fg=self.FG,
                              activebackground=self.BG_CARD,
                              activeforeground=self.ACCENT,
                              selectcolor=self.BG_SOFT,
                              relief="flat", bd=0,
                              font=("Segoe UI", 10), cursor="hand2",
                              highlightthickness=0)

    def _collect_settings_widgets(self, w):
        for c in w.winfo_children():
            if isinstance(c, (tk.Checkbutton, tk.Entry)):
                self._settings_widgets.append(c)
            self._collect_settings_widgets(c)

    def _set_settings_state(self, state):
        for c in self._settings_widgets:
            try:
                c.configure(state=state)
            except Exception:
                pass

    def start_work_panel(self, title):
        if self._work_panel_active:
            for w in self.right.winfo_children():
                w.destroy()
        else:
            self.right.pack(side="right", fill="y", padx=(14, 0))
            self._work_panel_active = True

        header = tk.Label(self.right, text=title, anchor="w",
                          bg=self.BG_CARD, fg=self.FG,
                          padx=20, pady=14,
                          font=("Segoe UI", 12, "bold"))
        header.pack(fill="x", side="top")

        sep = tk.Frame(self.right, height=1, bg=self.BORDER)
        sep.pack(fill="x", side="top")

        self._footer = tk.Frame(self.right, bg=self.BG_CARD)
        self._footer.pack(fill="x", side="bottom", pady=(6, 12), padx=14)

        self._footer_sep = tk.Frame(self.right, height=1, bg=self.BORDER)
        self._footer_sep.pack(fill="x", side="bottom")

        self._footer_help = tk.Label(self._footer, text="",
                                     justify="left", anchor="w",
                                     fg=self.FG_DIM, bg=self.BG_CARD,
                                     font=("Segoe UI", 9), wraplength=300)
        self._footer_help.pack(side="left", fill="x", expand=True)

        self._footer_btns = tk.Frame(self._footer, bg=self.BG_CARD)
        self._footer_btns.pack(side="right")

        body = tk.Frame(self.right, bg=self.BG_CARD, padx=14, pady=12)
        body.pack(fill="both", expand=True, side="top")

        return body

    def work_panel_footer(self, parent, on_ok=None, ok_text="OK",
                          cancel=False, cancel_value=None,
                          extra=None, help_text=None, cancel_text=None):
        if help_text:
            self._footer_help.config(text=help_text)

        if cancel_text is None:
            cancel_text = self._t("btn.cancel")

        def _cancel():
            self.stop_work_panel()
            self.bridge._ok(cancel_value)

        if extra:
            for text, cb in extra:
                RoundedButton(self._footer_btns, text, command=cb,
                              width=120, height=34, radius=17,
                              bg=self.BG_SOFT, fg=self.FG,
                              hover_bg=self.BG_HOVER, hover_fg=self.FG,
                              active_bg=self.BORDER,
                              font=("Segoe UI", 9)).pack(side="left", padx=4)

        if cancel:
            RoundedButton(self._footer_btns, cancel_text, command=_cancel,
                          width=110, height=34, radius=17,
                          bg="#3b2a2e", fg=self.ERR,
                          hover_bg="#4d3238", hover_fg=self.ERR,
                          active_bg="#5a3840",
                          font=("Segoe UI", 9)).pack(side="right", padx=4)

        if on_ok is not None and ok_text:
            RoundedButton(self._footer_btns, ok_text, command=on_ok,
                          width=150, height=34, radius=17,
                          bg=self.ACCENT, fg="#1a1c22",
                          hover_bg=self.ACCENT_D, hover_fg="#ffffff",
                          active_bg="#4a6fbf",
                          font=("Segoe UI", 9, "bold")).pack(side="right", padx=4)

    def stop_work_panel(self):
        if not self._work_panel_active:
            return
        for w in self.right.winfo_children():
            w.destroy()
        self.right.pack_forget()
        self._work_panel_active = False

    def append_log(self, text):
        t = str(text)
        tag = None
        low = t.lower()

        if ("ошибка" in low or "error" in low
                or "traceback" in low):
            tag = "err"
        elif ("[warn]" in low or "предупрежд" in low
              or "warning" in low or "⚠" in t):
            tag = "warn"
        elif ("анализ завершён" in low
              or "analysis complete" in low
              or "готово:" in low or "done:" in low):
            tag = "ok"
        elif (t.startswith("Обработка файла")
              or t.startswith("Processing file")
              or t.startswith("Итог")
              or t.startswith("File summary")
              or t.startswith("Разметка")
              or t.startswith("Soma marking")
              or t.startswith("Измерение")
              or t.startswith("Measuring")
              or t.startswith("Ручная")
              or t.startswith("Manual")
              or t.startswith("Средние")
              or t.startswith("Session averages")
              or t.startswith("Сегментация")
              or t.startswith("Segmenting")):
            tag = "accent"

        if tag:
            self.log_text.insert("end", t + "\n", tag)
        else:
            self.log_text.insert("end", t + "\n")
        self.log_text.see("end")
        self.root.update_idletasks()

    def _choose_file(self):
        if self.var_load_mode.get() == "one":
            p = filedialog.askopenfilename(
                filetypes=[
                    (self._t("dlg.files.images"),
                     "*.tif *.tiff *.nd2 *.czi *.lif *.ome.tif *.ome.tiff *.oib *.oif *.oir"),
                    ("TIFF", "*.tif *.tiff *.ome.tif *.ome.tiff"),
                    ("ND2", "*.nd2"),
                    ("CZI", "*.czi"),
                    ("LIF", "*.lif"),
                    ("Olympus", "*.oib *.oif *.oir"),
                    (self._t("dlg.files.all"), "*.*"),
                ])
            if p:
                self.var_file_path.set(p)
        else:
            d = filedialog.askdirectory()
            if d:
                self.var_file_path.set(d)

    def _read_metadata(self):
        self.append_log(self._t("log.metadata_reading", name="..."))
        try:
            path = self.var_file_path.get().strip()
            if not path or not os.path.exists(path):
                StyledDialog.show(self.root,
                                  self._t("ui.settings.metadata"),
                                  self._t("log.file_read_failed"),
                                  kind="warn")
                return

            if os.path.isdir(path):
                exts = ['.tif', '.tiff', '.nd2', '.ome.tif', '.ome.tiff',
                        '.oib', '.oif', '.oir']
                files = sorted(
                    os.path.join(path, f) for f in os.listdir(path)
                    if any(f.lower().endswith(e) for e in exts))
                if not files:
                    StyledDialog.show(self.root,
                                      self._t("ui.settings.metadata"),
                                      self._t("dlg.files.nothing"),
                                      kind="warn")
                    return
                sample = files[0]
            else:
                sample = path

            self.append_log(self._t("log.metadata_reading",
                                    name=os.path.basename(sample)))

            ctx = AnalysisContext()
            ctx.bridge = self.bridge
            ctx.lang = self.var_lang.get()
            file_io = FileIO(ctx)
            optics = Optics(ctx, file_io)

            try:
                meta = file_io._extract_voxel_size_from_file(sample)
            except Exception as e:
                self.append_log(self._t("log.metadata_scale_read_error", err=e))
                meta = {'xy_um': None, 'z_um': None, 'source': 'none'}

            if meta.get('xy_um'):
                self.var_xy.set(f"{meta['xy_um']:.4f}")
                self.lbl_xy_src.config(text=self._t("ui.settings.from_metadata"),
                                       fg=self.OK)
                self.append_log(self._t("log.metadata_xy_ok", v=meta['xy_um']))
            else:
                self.lbl_xy_src.config(text=self._t("ui.settings.default"),
                                       fg=self.FG_DIM)
                self.append_log(self._t("log.metadata_xy_none"))

            if meta.get('z_um'):
                self.var_z.set(f"{meta['z_um']:.4f}")
                self.lbl_z_src.config(text=self._t("ui.settings.from_metadata"),
                                      fg=self.OK)
                self.append_log(self._t("log.metadata_z_ok", v=meta['z_um']))
            else:
                self.lbl_z_src.config(text=self._t("ui.settings.default"),
                                      fg=self.FG_DIM)
                self.append_log(self._t("log.metadata_z_none"))

            try:
                op = optics._extract_optics_from_file(sample)
                channels = op.get('channels', [])
            except Exception as e:
                self.append_log(self._t("log.metadata_optics_read_error", err=e))
                channels = []

            cells_found = False
            dapi_found = False
            if channels:
                dapi_idx = op.get('dapi_channel_index')
                astro_idx = op.get('astro_channel_index')
                if astro_idx is not None and astro_idx < len(channels):
                    ch = channels[astro_idx]
                    if ch.get('na') is not None:
                        self.var_na_cells.set(f"{ch['na']:.3f}")
                    if ch.get('emission_um') is not None:
                        self.var_wl_cells.set(f"{ch['emission_um']*1000:.1f}")
                    if ch.get('na') is not None or ch.get('emission_um') is not None:
                        cells_found = True
                if dapi_idx is not None and dapi_idx < len(channels):
                    ch = channels[dapi_idx]
                    if ch.get('na') is not None:
                        self.var_na_dapi.set(f"{ch['na']:.3f}")
                    if ch.get('emission_um') is not None:
                        self.var_wl_dapi.set(f"{ch['emission_um']*1000:.1f}")
                    if ch.get('na') is not None or ch.get('emission_um') is not None:
                        dapi_found = True

            self.lbl_cells_src.config(
                text=(self._t("ui.settings.from_metadata") if cells_found
                      else self._t("ui.settings.default")),
                fg=self.OK if cells_found else self.FG_DIM)
            self.lbl_dapi_src.config(
                text=(self._t("ui.settings.from_metadata") if dapi_found
                      else self._t("ui.settings.default")),
                fg=self.OK if dapi_found else self.FG_DIM)

            self.append_log(self._t("log.metadata_fields_updated"))

        except Exception as e:
            self.append_log(self._t("log.error_prefix", err=e))
            self.append_log(traceback.format_exc())
            try:
                StyledDialog.show(self.root,
                                  self._t("dlg.error.title"),
                                  str(e),
                                  kind="error")
            except Exception:
                pass

    def start(self):
        if self.worker and self.worker.is_alive():
            StyledDialog.show(self.root,
                              self._t("app.title"),
                              self._t("log.analysis_done"),
                              kind="info")
            return
        path = self.var_file_path.get().strip()
        if not path or not os.path.exists(path):
            StyledDialog.show(self.root,
                              self._t("ui.settings.file"),
                              self._t("dlg.files.nothing"),
                              kind="warn")
            return
        self._switch_to_analysis_mode()
        self.btn_run.config_state("disabled")
        self._set_settings_state("disabled")
        self.worker = threading.Thread(target=self._run_core, daemon=True)
        self.worker.start()

    def _collect_settings(self):
        def f(s):
            return float(str(s).replace(",", ".").strip())
        return {
            "path": self.var_file_path.get().strip(),
            "load_mode": self.var_load_mode.get(),
            "cell_type": self.var_cell_type.get(),
            "image_mode": self.var_image_mode.get(),
            "use_deconv": bool(self.var_use_deconv.get()),
            "xy_um": f(self.var_xy.get()),
            "z_um":  f(self.var_z.get()),
            "na_cells": f(self.var_na_cells.get()),
            "wl_cells": f(self.var_wl_cells.get()),
            "eps_cells": f(self.var_eps_cells.get()),
            "na_dapi": f(self.var_na_dapi.get()),
            "wl_dapi": f(self.var_wl_dapi.get()),
            "eps_dapi": f(self.var_eps_dapi.get()),
            "do_filter": bool(self.var_do_filter.get()),
        }

    def _resolve_files(self, cfg):
        path = cfg["path"]
        exts = ['.tif', '.tiff', '.nd2', '.czi', '.lif',
                '.ome.tif', '.ome.tiff', '.oib', '.oif', '.oir']
        if os.path.isdir(path):
            return sorted(
                os.path.join(path, f) for f in os.listdir(path)
                if any(f.lower().endswith(e) for e in exts))
        return [path]

    def _run_core(self):
        try:
            try:
                cfg = self._collect_settings()
            except ValueError as e:
                self.root.after(0, lambda: StyledDialog.show(
                    self.root, self._t("dlg.error.title"),
                    str(e), kind="error"))
                self.root.after(0, lambda: self.btn_run.config_state("normal"))
                self.root.after(0, lambda: self._set_settings_state("normal"))
                return

            files = self._resolve_files(cfg)
            if not files:
                self.root.after(0, lambda: StyledDialog.show(
                    self.root, self._t("dlg.files.title_short"),
                    self._t("dlg.files.no_supported"), kind="warn"))
                return

            analyzer = Analyzer()
            analyzer.ctx.bridge = self.bridge
            analyzer.ctx.auto_mode = True
            analyzer.ctx.lang = self.var_lang.get()

            analyzer.ctx.cell_type = cfg["cell_type"]
            if cfg["cell_type"] == "astrocytes":
                analyzer.ctx.cell_type_ru = self._t("cell.astrocytes")
            elif cfg["cell_type"] == "microglia":
                analyzer.ctx.cell_type_ru = self._t("cell.microglia")
            else:
                analyzer.ctx.cell_type_ru = self._t("cell.neurons")

            analyzer.ctx.image_mode = cfg["image_mode"]
            analyzer.ctx.deconvolve_images = cfg["use_deconv"]

            analyzer.ctx.xy_scale_um = cfg["xy_um"]
            analyzer.ctx.z_scale_um = cfg["z_um"]
            analyzer.ctx.scale_initialized = True
            analyzer.ctx._input_xy_um = cfg["xy_um"]
            analyzer.ctx._input_z_um = cfg["z_um"]

            analyzer.ctx._deconv_params['na'] = cfg["na_cells"]
            analyzer.ctx._deconv_params['emission_wavelength_um'] = \
                cfg["wl_cells"] / 1000.0
            analyzer.ctx._deconv_params['wiener_epsilon'] = cfg["eps_cells"]
            analyzer.ctx._dapi_deconv_params['na'] = cfg["na_dapi"]
            analyzer.ctx._dapi_deconv_params['emission_wavelength_um'] = \
                cfg["wl_dapi"] / 1000.0
            analyzer.ctx._dapi_deconv_params['wiener_epsilon'] = cfg["eps_dapi"]

            analyzer.ctx._gui_channel_mode = "ask"
            analyzer.ctx._gui_do_filter = cfg["do_filter"]

            try:
                probe = FileIO(analyzer.ctx)._read_image_universal(files[0])
                if probe is not None:
                    if probe.ndim == 2:
                        y_dim, x_dim = probe.shape
                    elif probe.ndim in (3, 4):
                        y_dim, x_dim = probe.shape[-2:]
                    else:
                        y_dim, x_dim = 1024, 1024
                    analyzer.ctx.downsample = compute_downsample(y_dim, x_dim,
                                                                 target=512)
                    self.append_log(self._t("log.image_size",
                                            w=x_dim, h=y_dim,
                                            n=analyzer.ctx.downsample))
            except Exception as e:
                self.append_log(self._t("log.image_size_failed", err=e))

            self.append_log(self._t("log.scale_from_panel",
                                    xy=analyzer.ctx.xy_scale_um,
                                    z=analyzer.ctx.z_scale_um))

            analyzer.batch_process_files(files)
            self.append_log(self._t("log.analysis_done"))
        except Exception as e:
            self.append_log(self._t("log.error_prefix", err=e))
            self.append_log(traceback.format_exc())
        finally:
            self.root.after(0, lambda: self.btn_run.config_state("normal"))
            self.root.after(0, lambda: self._set_settings_state("normal"))

    def run(self):
        self.root.mainloop()

def main():
    App().run()


if __name__ == "__main__":
    main()