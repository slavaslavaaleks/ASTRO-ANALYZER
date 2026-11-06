# ASTRO-ANALYZER v1.0.2
# Copyright (C) 2026 Алексеев Вячеслав Иванович
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
    gaussian_filter, binary_dilation,
    median_filter, convolve, uniform_filter, sobel,
    distance_transform_edt,
    label, binary_fill_holes,
)

from skimage import measure
from skimage.transform import resize
from skimage.morphology import remove_small_objects

from skimage.filters import threshold_multiotsu
from scipy.ndimage import maximum_filter
from scipy.signal import find_peaks

warnings.filterwarnings('ignore')

TRANSLATIONS = {
    "ru": {
        "btn.yes": "Да",
        "btn.no": "Нет",

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
        "log.deconv_done": "    Деконволюция и ручной отбор завершены. Размер стека: {shape}",
        "log.deconv_dapi_done": "    DAPI-деконволюция завершена. Размер стека: {shape}",
        "log.soma_title": "\n    Разметка сом",
        "log.soma_no_points": "    Точки не добавлены. Сегментация отменена",
        "log.soma_points_set": "    Точек поставлено: {n}",
        "log.soma_cancelled": "    Расстановка отменена пользователем",
        "log.soma_starting": "    Сегментация сом...",
        "log.nuclei_title": "    Сегментация ядер...",
        "log.nuclei_no_dapi": "    DAPI-канал не используется — ядра не сегментируются",
        "log.nuclei_no_soma": "    Сомы не найдены",
        "log.nuclei_found": "    Ядер найдено: {n}",
        "log.cells_title": "    Сегментация клеток...",
        "log.cells_no_soma": "    Сомы не найдены",
        "log.cells_no_signal": "    После вычитания фона сигнал отсутствует",
        "log.cells_no_growth": "    Областей для роста не найдено",
        "log.cells_found": "    Клеток найдено: {n}",
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
        "log.viz_deconv_dapi_saved": "    Сохранено деконволюированное изображение DAPI: {shape}",
        "log.viz_mip_cells_failed": "    Не удалось сохранить MIP КЛЕТОК: {err}",
        "log.viz_mip_dapi_failed": "    Не удалось сохранить MIP DAPI: {err}",

        "log.granules_title": "    Сегментация гранул...",
        "log.viz_deconv_granules_saved": "    Сохранено деконволюированное изображение ГРАНУЛ: {shape}",
        "log.deconv_granules_title": "{name} — ГРАНУЛЫ (деконволюция, MIP)",
        "log.viz_mip_granules_failed": "    Не удалось сохранить MIP ГРАНУЛ: {err}",

        "log.slice_title": "Z-срез: {z}",
        "log.deconv_cells_title": "{name} — КЛЕТКИ (деконволюция, MIP)",
        "log.deconv_dapi_title": "{name} — DAPI (деконволюция, MIP)",
        "log.title_soma_slices": "Сегментация сом",
        "log.title_nuclei_slices": "Сегментация ядер",
        "log.title_cells_slices": "Визуализация целой клетки с отростками",
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
        "log.image_size": "    Размер изображения: {w}×{h} → рекомендован даунсемплинг в {n} раза",
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
        "log.filter_fused_correction": "    Пользователь указал на слипшиеся клетки. Возврат к коррекции сом.",
        "log.filter_artifacts_correction": "    Пользователь указал на артефакты. Возврат к удалению срезов.",
        "log.filter_artifacts_no_deconv": "    Пользователь указал на артефакты, но деконволюция была выключена. Требуется перезапуск.",
        "log.filter_cut_soft_mode": "    Пользователь указал на обрубленные отростки. Включен мягкий режим роста.",
        "log.filter_manual_deletion": "    Пользователь указал на лишние объекты. Запуск ручного удаления...",
        "log.soma_auto_starting": "    Автодетекция сом...",
        "log.soma_auto_found": "    Автоматически найдено сом: {n}",
        "log.soma_auto_none": "    Автодетекция не нашла сом — переход в ручной режим",
        "log.soma_auto_z_layers": "    Сом в (y={y}, x={x}): Z-слоёв {n}",

        "channel.cells": "клетки",
        "channel.dapi": "ядра (DAPI)",
        "channel.granules": "гранулы",
        "log.deconv_start_channel": "\n    Деконволюция канала: {ch}...",
        "log.deconv_done_channel": "    Деконволюция канала {ch} завершена. Размер стека: {shape}",
        "log.deconv_slices_deleted_channel": "    Канал {ch}: удалено срезов по вашему выбору: {n} (Z: {z})",
        "log.deconv_all_slices_kept_channel": "    Канал {ch}: все срезы оставлены",
        "dlg.slices.title_selector_channel": "Выберите срезы для удаления — {ch}.",

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
        "ui.settings.granules": "Гранулы:",
        "ui.settings.default": "[по умолчанию]",
        "ui.settings.from_metadata": "[из метаданных]",
        "ui.settings.filter": "Ручная фильтрация клеток:",
        "ui.settings.run_filter": "Запускать после сегментации",
        "ui.settings.run": "▶   Запустить анализ",
        "ui.settings.exit": "Выход",
        "ui.settings.back_to_log": "→ Вернуться к логу",
        "ui.settings.back_to_settings": "← Настройки",
        "ui.settings.log_title": "Протокол работы",
        "ui.settings.na": "NA:",
        "ui.settings.lambda_nm": "λ(нм):",
        "ui.settings.epsilon": "ε:",
        "ui.settings.auto_soma": "Автодетекция сом:",
        "ui.settings.auto_soma_check": "Определять сомы автоматически [Beta]",
        "ui.header.lang": "Язык:",

        "ui.settings.not_in_metadata": "нет в метаданных: {missing}",
        "ui.settings.meta_na": "NA",
        "ui.settings.meta_wl": "λ",

        "btn.ok": "OK",
        "btn.cancel": "Отмена",
        "btn.apply": "Применить",
        "btn.done": "Готово",
        "btn.continue": "Продолжить",
        "btn.skip": "Пропустить",
        "btn.start_analysis": "Запустить сегментацию",
        "btn.choose": "Выбрать…",
        "btn.delete_return_slice": "Удалить / Вернуть срез",
        "btn.to_manual_soma": "К ручной расстановке сом",

        "dlg.files.title": "Выбор файлов",
        "dlg.files.one_or_dir": "Один файл или папка?",
        "dlg.files.one": "Один файл",
        "dlg.files.dir": "Папка с файлами",
        "dlg.files.selected_n": "Выбрано: {n} файл(ов)",
        "dlg.files.nothing": "Ничего не выбрано.",
        "dlg.files.images": "Изображения",
        "dlg.files.all": "Все файлы",

        "dlg.preprocess.title": "Предобработка изображения",
        "dlg.preprocess.size_line": "Размер изображения: {w}×{h} px",
        "dlg.preprocess.ds_question": "Проводить даунсемплинг до размера 512 px по максимальному измерению?",
        "dlg.preprocess.crop_question": "Хотите обрезать изображение по Oxy?",
        "dlg.preprocess.help": "Отметьте нужные пункты и нажмите «Применить».",
        "dlg.crop.title": "Обрезка изображения (Oxy)",
        "dlg.crop.width_um": "Ширина (мкм):",
        "dlg.crop.height_um": "Высота (мкм):",
        "dlg.crop.coords": "X: [{x0:.1f}..{x1:.1f}]  Y: [{y0:.1f}..{y1:.1f}]   W={w:.1f}  H={h:.1f} (мкм)",
        "dlg.crop.help": "Перетащите рамку или её ручки. Размер можно ввести в поля и нажать Enter.",

        "dlg.channels.title": "Выбор каналов",
        "dlg.channels.instructions": "Отметьте каналы, которые войдут в анализ",
        "dlg.channels.channel_n": "КАНАЛ {n}",
        "dlg.channels.choose_one": "Выберите хотя бы один канал.",

        "dlg.channels.role_cells": "Клетки",
        "dlg.channels.role_dapi": "DAPI",
        "dlg.channels.role_granules": "Гранулы",
        "dlg.channels.role_none": "Не включать в анализ",
        "dlg.channels.role_used": "Уже назначен другому каналу",
        "dlg.channels.help": (
            "Для каждого канала выберите роль из списка под картинкой.\n"
            "«Клетки» — обязательная роль.\n"
            "«DAPI» — ядра (необязательно).\n"
            "«Гранулы» — внутриклеточные гранулы (необязательно).\n"
            "«Не включать в анализ» — канал игнорируется.\n"
            "Роль, уже занятая другим каналом, из списка исчезает."
        ),

        "log.granule_channel": "    Гранулы: канал {idx}",

        "dlg.metrics.title": "Метрики для измерения",
        "dlg.metrics.group_targets": "Что измерять",
        "dlg.metrics.group_shape_2d_cell": "Форма клетки (2D)",
        "dlg.metrics.group_shape_2d_soma": "Форма сомы (2D)",
        "dlg.metrics.group_shape_2d_nucleus": "Форма ядра (2D)",
        "dlg.metrics.group_shape_3d_cell": "Форма клетки (3D)",
        "dlg.metrics.group_shape_3d_soma": "Форма сомы (3D)",
        "dlg.metrics.group_shape_3d_nucleus": "Форма ядра (3D)",
        "dlg.metrics.group_intensity_cell": "Яркость клетки",
        "dlg.metrics.group_intensity_soma": "Яркость сомы",
        "dlg.metrics.group_intensity_nucleus": "Яркость ядра",
        "dlg.metrics.group_texture_cell": "Текстура клетки",
        "dlg.metrics.group_texture_soma": "Текстура сомы",
        "dlg.metrics.group_texture_nucleus": "Текстура ядра",
        "dlg.metrics.group_nucleus": "Ядро и цитоплазма",
        "dlg.metrics.group_skeleton": "Скелет / отростки",
        "dlg.metrics.group_granules": "Гранулы",
        "dlg.metrics.select_all": "Выбрать все",
        "dlg.metrics.select_none": "Снять все",
        "dlg.metrics.btn": "Метрики…",
        "dlg.metrics.hint": "Отметьте метрики для CSV. Каждая группа применяется только к своей структуре.",

        "dlg.metrics.item.target_cell": "Клетка",
        "dlg.metrics.item.target_soma": "Сома",
        "dlg.metrics.item.target_nucleus": "Ядро",

        "dlg.metrics.item.cell_area": "Площадь клетки",
        "dlg.metrics.item.cell_perimeter": "Периметр клетки",
        "dlg.metrics.item.cell_solidity": "Солидность клетки",
        "dlg.metrics.item.cell_elongation": "Вытянутость клетки",
        "dlg.metrics.item.cell_form_factor": "Фактор формы клетки",
        "dlg.metrics.item.cell_roundness": "Округлость клетки",
        "dlg.metrics.item.cell_aspect_ratio": "Соотношение сторон клетки",
        "dlg.metrics.item.cell_eccentricity": "Эксцентриситет клетки",
        "dlg.metrics.item.cell_convexity": "Выпуклость клетки",
        "dlg.metrics.item.cell_extent": "Заполнение клетки",
        "dlg.metrics.item.cell_feret": "Диаметры Фере клетки",
        "dlg.metrics.item.cell_major_minor": "Оси клетки",
        "dlg.metrics.item.cell_bbox_hull": "Габариты и оболочка клетки",
        "dlg.metrics.item.cell_moments": "Моменты формы клетки",
        "dlg.metrics.item.cell_equiv_diam": "Эквивалентный диаметр клетки",

        "dlg.metrics.item.soma_area": "Площадь сомы",
        "dlg.metrics.item.soma_perimeter": "Периметр сомы",
        "dlg.metrics.item.soma_solidity": "Солидность сомы",
        "dlg.metrics.item.soma_elongation": "Вытянутость сомы",
        "dlg.metrics.item.soma_form_factor": "Фактор формы сомы",
        "dlg.metrics.item.soma_roundness": "Округлость сомы",
        "dlg.metrics.item.soma_eccentricity": "Эксцентриситет сомы",
        "dlg.metrics.item.soma_convexity": "Выпуклость сомы",

        "dlg.metrics.item.nucleus_area": "Площадь ядра",
        "dlg.metrics.item.nucleus_perimeter": "Периметр ядра",
        "dlg.metrics.item.nucleus_solidity": "Солидность ядра",
        "dlg.metrics.item.nucleus_elongation": "Вытянутость ядра",
        "dlg.metrics.item.nucleus_form_factor": "Фактор формы ядра",
        "dlg.metrics.item.nucleus_roundness": "Округлость ядра",
        "dlg.metrics.item.nucleus_eccentricity": "Эксцентриситет ядра",
        "dlg.metrics.item.nucleus_convexity": "Выпуклость ядра",

        "dlg.metrics.item.cell_volume": "Объём клетки",
        "dlg.metrics.item.cell_surface": "Площадь поверхности клетки",
        "dlg.metrics.item.cell_sphericity": "Сферичность клетки",
        "dlg.metrics.item.cell_sav": "Площадь/Объём клетки",
        "dlg.metrics.item.cell_convexity_3d": "Выпуклость 3D клетки",
        "dlg.metrics.item.cell_extent_3d": "Заполнение 3D клетки",
        "dlg.metrics.item.cell_axes_3d": "Главные оси 3D клетки",

        "dlg.metrics.item.soma_volume": "Объём сомы",
        "dlg.metrics.item.soma_elongation_3d": "Вытянутость 3D сомы",
        "dlg.metrics.item.soma_flatness": "Плоскостность сомы",
        "dlg.metrics.item.soma_extent_3d": "Заполнение 3D сомы",
        "dlg.metrics.item.soma_convexity_3d": "Выпуклость 3D сомы",

        "dlg.metrics.item.nucleus_volume": "Объём ядра",
        "dlg.metrics.item.nucleus_elongation_3d": "Вытянутость 3D ядра",
        "dlg.metrics.item.nucleus_flatness": "Плоскостность ядра",
        "dlg.metrics.item.nucleus_extent_3d": "Заполнение 3D ядра",
        "dlg.metrics.item.nucleus_convexity_3d": "Выпуклость 3D ядра",

        "dlg.metrics.item.cell_mean_intensity": "Средняя яркость клетки",
        "dlg.metrics.item.cell_intensity_stats": "Разброс яркости клетки (Std/CV/IQR)",
        "dlg.metrics.item.cell_intensity_percentiles": "Перцентили яркости клетки",
        "dlg.metrics.item.cell_intensity_entropy_gini": "Энтропия и коэффициент Джини клетки",
        "dlg.metrics.item.cell_intensity_highlow": "Доля ярких/тёмных пикселей клетки",
        "dlg.metrics.item.cell_snr_cnr_sbr": "Соотношения сигнал/шум (SNR/CNR/SBR)",

        "dlg.metrics.item.soma_mean_intensity": "Средняя яркость сомы",
        "dlg.metrics.item.soma_intensity_stats": "Разброс яркости сомы (Std/CV)",
        "dlg.metrics.item.soma_intensity_entropy": "Энтропия яркости сомы",

        "dlg.metrics.item.nucleus_mean_intensity": "Средняя яркость ядра",
        "dlg.metrics.item.nucleus_intensity_stats": "Разброс яркости ядра (Std/CV)",
        "dlg.metrics.item.nucleus_intensity_entropy": "Энтропия яркости ядра",

        "dlg.metrics.item.cell_local_variance": "Локальная дисперсия клетки",
        "dlg.metrics.item.cell_glcm": "Текстурные признаки клетки (GLCM)",
        "dlg.metrics.item.cell_fractal": "Фрактальная размерность клетки",
        "dlg.metrics.item.soma_glcm": "Текстурные признаки сомы (GLCM)",
        "dlg.metrics.item.nucleus_glcm": "Текстурные признаки ядра (GLCM)",

        "dlg.metrics.item.cell_nuc_cyto": "Отношение ядро/цитоплазма",
        "dlg.metrics.item.cell_nuc_intensity": "Соотношение яркости ядра и цитоплазмы",
        "dlg.metrics.item.cell_centroid_dist": "Смещение центроидов (ядро и сома)",
        "dlg.metrics.item.cell_process": "Количество отростков",
        "dlg.metrics.item.cell_skeleton": "Скелет: длина, ветвления, извитость",

        "dlg.metrics.item.granule_intensity": "Яркость гранул",
        "dlg.metrics.item.granule_size_count": "Размер и количество гранул",
        "dlg.metrics.item.granule_shape": "Форма гранул",
        "dlg.metrics.item.granule_spatial": "Пространственное распределение гранул",
        "dlg.metrics.item.granule_sbr": "Контраст гранул (SBR)",

        "dlg.metrics.group_targets_legacy": "К каким структурам применять форму/текстуру",

        "dlg.dapi.title": "Определение DAPI-канала",
        "dlg.dapi.text": "Автоматически определить не удалось.\nКликните по каналу, который является DAPI:",
        "dlg.dapi.help": "Кликните по каналу с ядрами.",

        "dlg.soma.title": "Расстановка сом",
        "dlg.soma.counter": "Точек: {n}",
        "dlg.soma.coords": "XY = ({y}, {x})",
        "dlg.soma.help": "ЛКМ — добавить сому,\nCtrl+ЛКМ — удалить ближайшую,\nПКМ — открыть лупу (оверлей).\nВ оверлее: слайдер Z, ЛКМ — поставить сому.",
        "dlg.soma.overlay_z": "Z = {z}",

        "dlg.soma.type_label": "Тип новой сомы:",
        "dlg.soma.type_solid": "Обычная",
        "dlg.soma.type_donut": "Бублик",
        "dlg.soma.help": "ЛКМ — добавить сому (тип задается переключателем),\nShift+ЛКМ — изменить тип ближайшей сомы,\nCtrl+ЛКМ — удалить ближайшую,\nПКМ — открыть лупу (оверлей).",

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
        "dlg.slices.title_selector": "Выберите срезы для удаления.",

        "dlg.filter.prob_not_cells": "Нужно удалить часть сегментов",
        "dlg.filter.warn_delete": "Сейчас откроется экран удаления. Кликните по объектам, которые нужно удалить из анализа.",
        "dlg.filter.delete_title": "Удаление лишних сегментов",
        "dlg.filter.delete_help": "Кликните по объекту, чтобы удалить его из анализа.",
        "dlg.filter.quality_title": "Качество сегментации",
        "dlg.filter.quality_ask": "Устраивает ли вас качество текущей сегментации клеток?",
        "dlg.filter.problem_title": "В чём проблема?",
        "dlg.filter.problem_ask": "Укажите, что пошло не так (это поможет перенаправить вас в нужное место):",
        "dlg.filter.prob_fused": "Нужно поправить сомы",
        "dlg.filter.prob_artifacts": "Нужно убрать лишние Z-срезы",
        "dlg.filter.prob_cut": "Нужно увеличить агрессивность роста",
        "dlg.filter.warn_soma": "Сейчас откроется экран коррекции сом. Добавьте пропущенные сомы или удалите лишние, чтобы разделить слипшиеся клетки.",
        "dlg.filter.warn_slices": "Сейчас откроется экран удаления Z-срезов. Уберите срезы, на которых видны артефакты или шум.",
        "dlg.filter.warn_slices_fallback": "Деконволюция была выключена, поэтому вам нужно вернуться в настройки и перезапустить анализ, чтобы убрать артефакты.",
        "dlg.filter.warn_cut": "Рост отростков будет усилен (снижен порог контраста). Если это не поможет, попробуйте в следующий раз выбрать другой тип изображения или уменьшить деконволюцию.",
        "dlg.filter.resegmenting_soft": "    Применяется мягкая коррекция роста отростков...",

        "entity.cells": "клеток",
        "entity.soma": "сом",
        "entity.nuclei": "ядер",
        "entity.granules": "гранул",

        "dlg.filter.quality_ask_entity": "Устраивает ли вас качество текущей сегментации {entity}?",
        "dlg.filter.prob_not_cells": "Нужно удалить часть сегментов",

        "dlg.unknown": "Неизвестный диалог: {kind}",
        "log.dialog_error": "ОШИБКА ДИАЛОГА: {err}",
        "log.donut_volume_shape_error": "Ожидается astro_norm_3d формы (Z, Y, X)",

        "dlg.error.title": "Ошибка",
        "dlg.files.title_short": "Файлы",
        "dlg.files.no_supported": "Поддерживаемые файлы не найдены.",

        "app.title": "ASTRO-ANALYZER",
        "app.header_title": "3D анализ клеток нервной ткани",
    },

    "en": {
        "btn.yes": "Yes",
        "btn.no": "No",
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
        "log.nuclei_title": "    Segmenting nuclei...",
        "log.nuclei_no_dapi": "    DAPI channel not used — nuclei are not segmented",
        "log.nuclei_no_soma": "    No somas found",
        "log.nuclei_found": "    Nuclei found: {n}",
        "log.cells_title": "    Segmenting cells...",
        "log.cells_no_soma": "    No somas found",
        "log.cells_no_signal": "    No signal left after background subtraction",
        "log.cells_no_growth": "    No regions available for growth",
        "log.cells_found": "    Cells found: {n}",
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
        "log.viz_deconv_dapi_saved": "    Deconvolved DAPI image saved: {shape}",
        "log.viz_mip_cells_failed": "    Failed to save CELLS MIP: {err}",
        "log.viz_mip_dapi_failed": "    Failed to save DAPI MIP: {err}",
        "log.slice_title": "Z-slice: {z}",

        "log.granules_title": "    Segmenting granules...",
        "log.viz_deconv_granules_saved": "    Deconvolved GRANULES image saved: {shape}",
        "log.deconv_granules_title": "{name} — GRANULES (deconvolution, MIP)",
        "log.viz_mip_granules_failed": "    Failed to save GRANULES MIP: {err}",

        "log.deconv_cells_title": "{name} — CELLS (deconvolution, MIP)",
        "log.deconv_dapi_title": "{name} — DAPI (deconvolution, MIP)",
        "log.title_soma_slices": "Soma segmentation",
        "log.title_nuclei_slices": "Nuclei segmentation",
        "log.title_cells_slices": "Whole cell visualization with processes",
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
        "log.image_size": "    Image size: {w}×{h} → downsampling by {n} is recommended",
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
        "log.filter_fused_correction": "    User reported merged cells. Returning to soma correction.",
        "log.filter_artifacts_correction": "    User reported artifacts. Returning to slice removal.",
        "log.filter_artifacts_no_deconv": "    User reported artifacts, but deconvolution was disabled. Restart required.",
        "log.filter_cut_soft_mode": "    User reported truncated processes. Soft growth mode enabled.",
        "log.filter_manual_deletion": "    User reported extra objects. Starting manual deletion...",
        "log.soma_auto_starting": "    Auto-detecting somas...",
        "log.soma_auto_found": "    Somas auto-detected: {n}",
        "log.soma_auto_none": "    Auto-detection found no somas — falling back to manual mode",
        "log.soma_auto_z_layers": "    Soma at (y={y}, x={x}): Z-layers {n}",

        "channel.cells": "cells",
        "channel.dapi": "nuclei (DAPI)",
        "channel.granules": "granules",
        "log.deconv_start_channel": "\n    Deconvolving channel: {ch}...",
        "log.deconv_done_channel": "    Channel {ch} deconvolution complete. Stack shape: {shape}",
        "log.deconv_slices_deleted_channel": "    Channel {ch}: slices removed by your choice: {n} (Z: {z})",
        "log.deconv_all_slices_kept_channel": "    Channel {ch}: all slices kept",
        "dlg.slices.title_selector_channel": "Choose slices to delete — {ch}.",

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
        "ui.settings.granules": "Granules:",
        "ui.settings.default": "[default]",
        "ui.settings.from_metadata": "[from metadata]",
        "ui.settings.filter": "Manual cell filtering:",
        "ui.settings.run_filter": "Run after segmentation",
        "ui.settings.run": "▶   Run analysis",
        "ui.settings.exit": "Exit",
        "ui.settings.back_to_log": "→ Back to log",
        "ui.settings.back_to_settings": "← Settings",
        "ui.settings.log_title": "Log",
        "ui.settings.na": "NA:",
        "ui.settings.lambda_nm": "λ(nm):",
        "ui.settings.epsilon": "ε:",
        "ui.settings.auto_soma": "Soma auto-detection:",
        "ui.settings.auto_soma_check": "Detect somas automatically [Beta]",
        "ui.header.lang": "Lang:",

        "ui.settings.not_in_metadata": "not in metadata: {missing}",
        "ui.settings.meta_na": "NA",
        "ui.settings.meta_wl": "λ",

        "btn.ok": "OK",
        "btn.cancel": "Cancel",
        "btn.apply": "Apply",
        "btn.done": "Done",
        "btn.continue": "Continue",
        "btn.skip": "Skip",
        "btn.start_analysis": "Start segmentation",
        "btn.choose": "Choose…",
        "btn.delete_return_slice": "Delete / Restore slice",
        "btn.to_manual_soma": "To manual soma marking",

        "dlg.files.title": "File selection",
        "dlg.files.one_or_dir": "Single file or folder?",
        "dlg.files.one": "Single file",
        "dlg.files.dir": "Folder with files",
        "dlg.files.selected_n": "Selected: {n} file(s)",
        "dlg.files.nothing": "Nothing selected.",
        "dlg.files.images": "Images",
        "dlg.files.all": "All files",

        "dlg.preprocess.title": "Image preprocessing",
        "dlg.preprocess.size_line": "Image size: {w}×{h} px",
        "dlg.preprocess.ds_question": "Downsample to 512 px along the largest dimension?",
        "dlg.preprocess.crop_question": "Do you want to crop the image in the Oxy plane?",
        "dlg.preprocess.help": "Tick the required options and press \"Apply\".",
        "dlg.crop.title": "Image crop (Oxy)",
        "dlg.crop.width_um": "Width (µm):",
        "dlg.crop.height_um": "Height (µm):",
        "dlg.crop.coords": "X: [{x0:.1f}..{x1:.1f}]  Y: [{y0:.1f}..{y1:.1f}]   W={w:.1f}  H={h:.1f} (µm)",
        "dlg.crop.help": "Drag the box or its handles. You can enter the size and press Enter.",

        "dlg.channels.title": "Channel selection",
        "dlg.channels.instructions": "Mark the channels to include in the analysis",
        "dlg.channels.channel_n": "CHANNEL {n}",
        "dlg.channels.choose_one": "Select at least one channel.",

        "dlg.channels.role_cells": "Cells",
        "dlg.channels.role_dapi": "DAPI",
        "dlg.channels.role_granules": "Granules",
        "dlg.channels.role_none": "Not in analysis",
        "dlg.channels.role_used": "Already assigned to another channel",
        "dlg.channels.help": (
            "Pick a role for each channel from the dropdown below the thumbnail.\n"
            "«Cells» — required role.\n"
            "«DAPI» — nuclei (optional).\n"
            "«Granules» — intracellular granules (optional).\n"
            "«Not in analysis» — channel ignored.\n"
            "A role already taken by another channel disappears from the list."
        ),

        "log.granule_channel": "    Granules: channel {idx}",

        "dlg.metrics.title": "Metrics to measure",
        "dlg.metrics.group_targets": "What to measure",
        "dlg.metrics.group_shape_2d_cell": "Cell shape (2D)",
        "dlg.metrics.group_shape_2d_soma": "Soma shape (2D)",
        "dlg.metrics.group_shape_2d_nucleus": "Nucleus shape (2D)",
        "dlg.metrics.group_shape_3d_cell": "Cell shape (3D)",
        "dlg.metrics.group_shape_3d_soma": "Soma shape (3D)",
        "dlg.metrics.group_shape_3d_nucleus": "Nucleus shape (3D)",
        "dlg.metrics.group_intensity_cell": "Cell intensity",
        "dlg.metrics.group_intensity_soma": "Soma intensity",
        "dlg.metrics.group_intensity_nucleus": "Nucleus intensity",
        "dlg.metrics.group_texture_cell": "Cell texture",
        "dlg.metrics.group_texture_soma": "Soma texture",
        "dlg.metrics.group_texture_nucleus": "Nucleus texture",
        "dlg.metrics.group_nucleus": "Nucleus and cytoplasm",
        "dlg.metrics.group_skeleton": "Skeleton / processes",
        "dlg.metrics.group_granules": "Granules",
        "dlg.metrics.select_all": "Select all",
        "dlg.metrics.select_none": "Deselect all",
        "dlg.metrics.btn": "Metrics…",
        "dlg.metrics.hint": "Tick metrics to include in CSV. Each group applies only to its own structure.",

        "dlg.metrics.item.target_cell": "Cell",
        "dlg.metrics.item.target_soma": "Soma",
        "dlg.metrics.item.target_nucleus": "Nucleus",

        "dlg.metrics.item.cell_area": "Cell area",
        "dlg.metrics.item.cell_perimeter": "Cell perimeter",
        "dlg.metrics.item.cell_solidity": "Cell solidity",
        "dlg.metrics.item.cell_elongation": "Cell elongation",
        "dlg.metrics.item.cell_form_factor": "Cell form factor",
        "dlg.metrics.item.cell_roundness": "Cell roundness",
        "dlg.metrics.item.cell_aspect_ratio": "Cell aspect ratio",
        "dlg.metrics.item.cell_eccentricity": "Cell eccentricity",
        "dlg.metrics.item.cell_convexity": "Cell convexity",
        "dlg.metrics.item.cell_extent": "Cell extent",
        "dlg.metrics.item.cell_feret": "Cell Feret diameters",
        "dlg.metrics.item.cell_major_minor": "Cell major/minor axes",
        "dlg.metrics.item.cell_bbox_hull": "Cell bbox and convex hull",
        "dlg.metrics.item.cell_moments": "Cell shape moments",
        "dlg.metrics.item.cell_equiv_diam": "Cell equivalent diameter",

        "dlg.metrics.item.soma_area": "Soma area",
        "dlg.metrics.item.soma_perimeter": "Soma perimeter",
        "dlg.metrics.item.soma_solidity": "Soma solidity",
        "dlg.metrics.item.soma_elongation": "Soma elongation",
        "dlg.metrics.item.soma_form_factor": "Soma form factor",
        "dlg.metrics.item.soma_roundness": "Soma roundness",
        "dlg.metrics.item.soma_eccentricity": "Soma eccentricity",
        "dlg.metrics.item.soma_convexity": "Soma convexity",

        "dlg.metrics.item.nucleus_area": "Nucleus area",
        "dlg.metrics.item.nucleus_perimeter": "Nucleus perimeter",
        "dlg.metrics.item.nucleus_solidity": "Nucleus solidity",
        "dlg.metrics.item.nucleus_elongation": "Nucleus elongation",
        "dlg.metrics.item.nucleus_form_factor": "Nucleus form factor",
        "dlg.metrics.item.nucleus_roundness": "Nucleus roundness",
        "dlg.metrics.item.nucleus_eccentricity": "Nucleus eccentricity",
        "dlg.metrics.item.nucleus_convexity": "Nucleus convexity",

        "dlg.metrics.item.cell_volume": "Cell volume",
        "dlg.metrics.item.cell_surface": "Cell surface area",
        "dlg.metrics.item.cell_sphericity": "Cell sphericity",
        "dlg.metrics.item.cell_sav": "Cell surface/volume",
        "dlg.metrics.item.cell_convexity_3d": "Cell 3D convexity",
        "dlg.metrics.item.cell_extent_3d": "Cell 3D extent",
        "dlg.metrics.item.cell_axes_3d": "Cell 3D principal axes",

        "dlg.metrics.item.soma_volume": "Soma volume",
        "dlg.metrics.item.soma_elongation_3d": "Soma 3D elongation",
        "dlg.metrics.item.soma_flatness": "Soma flatness",
        "dlg.metrics.item.soma_extent_3d": "Soma 3D extent",
        "dlg.metrics.item.soma_convexity_3d": "Soma 3D convexity",

        "dlg.metrics.item.nucleus_volume": "Nucleus volume",
        "dlg.metrics.item.nucleus_elongation_3d": "Nucleus 3D elongation",
        "dlg.metrics.item.nucleus_flatness": "Nucleus flatness",
        "dlg.metrics.item.nucleus_extent_3d": "Nucleus 3D extent",
        "dlg.metrics.item.nucleus_convexity_3d": "Nucleus 3D convexity",

        "dlg.metrics.item.cell_mean_intensity": "Cell mean intensity",
        "dlg.metrics.item.cell_intensity_stats": "Cell intensity spread (Std/CV/IQR)",
        "dlg.metrics.item.cell_intensity_percentiles": "Cell intensity percentiles",
        "dlg.metrics.item.cell_intensity_entropy_gini": "Cell entropy and Gini",
        "dlg.metrics.item.cell_intensity_highlow": "Cell high/low fraction",
        "dlg.metrics.item.cell_snr_cnr_sbr": "SNR/CNR/SBR",

        "dlg.metrics.item.soma_mean_intensity": "Soma mean intensity",
        "dlg.metrics.item.soma_intensity_stats": "Soma intensity spread (Std/CV)",
        "dlg.metrics.item.soma_intensity_entropy": "Soma intensity entropy",

        "dlg.metrics.item.nucleus_mean_intensity": "Nucleus mean intensity",
        "dlg.metrics.item.nucleus_intensity_stats": "Nucleus intensity spread (Std/CV)",
        "dlg.metrics.item.nucleus_intensity_entropy": "Nucleus intensity entropy",

        "dlg.metrics.item.cell_local_variance": "Cell local variance",
        "dlg.metrics.item.cell_glcm": "Cell GLCM texture",
        "dlg.metrics.item.cell_fractal": "Cell fractal dimension",
        "dlg.metrics.item.soma_glcm": "Soma GLCM texture",
        "dlg.metrics.item.nucleus_glcm": "Nucleus GLCM texture",

        "dlg.metrics.item.cell_nuc_cyto": "Nucleus/cytoplasm ratio",
        "dlg.metrics.item.cell_nuc_intensity": "Nucleus/cytoplasm intensity ratio",
        "dlg.metrics.item.cell_centroid_dist": "Centroid offset (nucleus and soma)",
        "dlg.metrics.item.cell_process": "Process count",
        "dlg.metrics.item.cell_skeleton": "Skeleton: length, branches, tortuosity",

        "dlg.metrics.item.granule_intensity": "Granule intensity",
        "dlg.metrics.item.granule_size_count": "Granule size and count",
        "dlg.metrics.item.granule_shape": "Granule shape",
        "dlg.metrics.item.granule_spatial": "Granule spatial distribution",
        "dlg.metrics.item.granule_sbr": "Granule SBR",

        "dlg.metrics.group_targets_legacy": "Apply shape/texture to",

        "dlg.dapi.title": "DAPI channel detection",
        "dlg.dapi.text": "Automatic detection failed.\nClick the channel that contains DAPI:",
        "dlg.dapi.help": "Click the channel with nuclei.",

        "dlg.soma.title": "Soma marking",
        "dlg.soma.counter": "Points: {n}",
        "dlg.soma.coords": "XY = ({y}, {x})",
        "dlg.soma.help": "LMB — add soma,\nCtrl+LMB — delete nearest,\nRMB — open magnifier (overlay).\nIn overlay: Z slider, LMB — place soma.",
        "dlg.soma.overlay_z": "Z = {z}",

        "dlg.soma.type_label": "New soma type:",
        "dlg.soma.type_solid": "Solid",
        "dlg.soma.type_donut": "Donut",
        "dlg.soma.help": "LMB — add soma (type chosen below),\nShift+LMB — toggle type of nearest,\nCtrl+LMB — delete nearest,\nRMB — open magnifier (overlay).",

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
        "dlg.slices.title_selector": "Choose slices to delete.",

        "dlg.filter.prob_not_cells": "Need to delete some segments",
        "dlg.filter.warn_delete": "Deletion screen will open. Click on the objects you want to remove from the analysis.",
        "dlg.filter.delete_title": "Remove extra segments",
        "dlg.filter.delete_help": "Click on an object to remove it from the analysis.",
        "dlg.filter.quality_title": "Segmentation quality",
        "dlg.filter.quality_ask": "Are you satisfied with the current cell segmentation quality?",
        "dlg.filter.problem_title": "What is the problem?",
        "dlg.filter.problem_ask": "Please specify what went wrong (this will redirect you to the right step):",
        "dlg.filter.prob_fused": "Need to fix somas",
        "dlg.filter.prob_artifacts": "Need to remove extra Z-slices",
        "dlg.filter.prob_cut": "Need to increase growth aggressiveness",
        "dlg.filter.warn_soma": "Soma correction screen will open. Add missed somas or delete extra ones to split merged cells.",
        "dlg.filter.warn_slices": "Z-slice selection screen will open. Remove slices with artifacts or noise.",
        "dlg.filter.warn_slices_fallback": "Deconvolution was disabled, so you need to go back to settings and restart the analysis to remove artifacts.",
        "dlg.filter.warn_cut": "Process growth will be boosted (contrast threshold lowered). If it doesn't help, try another image type or reduce deconvolution next time.",
        "dlg.filter.resegmenting_soft": "    Applying soft process growth correction...",

        "entity.cells": "cells",
        "entity.soma": "somas",
        "entity.nuclei": "nuclei",
        "entity.granules": "granules",
        "dlg.filter.quality_ask_entity": "Are you satisfied with the current segmentation quality of {entity}?",
        "dlg.filter.prob_not_cells": "Need to delete some segments",

        "dlg.unknown": "Unknown dialog: {kind}",
        "log.dialog_error": "DIALOG ERROR: {err}",
        "log.donut_volume_shape_error": "Expected astro_norm_3d of shape (Z, Y, X)",

        "dlg.error.title": "Error",
        "dlg.files.title_short": "Files",
        "dlg.files.no_supported": "No supported files found.",

        "app.title": "ASTRO-ANALYZER",
        "app.header_title": "3D analysis of nerve tissue cells",
    }
}

_2D_ONLY_METRICS = {
    "cell_area", "cell_perimeter", "cell_solidity", "cell_elongation",
    "cell_form_factor", "cell_roundness", "cell_aspect_ratio",
    "cell_eccentricity", "cell_convexity", "cell_extent",
    "cell_feret", "cell_major_minor", "cell_bbox_hull",
    "cell_moments", "cell_equiv_diam",
    "soma_area", "soma_perimeter", "soma_solidity", "soma_elongation",
    "soma_form_factor", "soma_roundness", "soma_eccentricity",
    "soma_convexity",
    "nucleus_area", "nucleus_perimeter", "nucleus_solidity",
    "nucleus_elongation", "nucleus_form_factor", "nucleus_roundness",
    "nucleus_eccentricity", "nucleus_convexity",
}

_3D_ONLY_METRICS = {
    "cell_volume", "cell_surface", "cell_sphericity", "cell_sav",
    "cell_convexity_3d", "cell_extent_3d", "cell_axes_3d",
    "soma_volume", "soma_elongation_3d", "soma_flatness",
    "soma_extent_3d", "soma_convexity_3d",
    "nucleus_volume", "nucleus_elongation_3d", "nucleus_flatness",
    "nucleus_extent_3d", "nucleus_convexity_3d",
}

METRIC_GROUPS = [
    ("group_shape_2d_cell", [
        "cell_area", "cell_perimeter", "cell_solidity", "cell_elongation",
        "cell_form_factor", "cell_roundness", "cell_aspect_ratio",
        "cell_eccentricity", "cell_convexity", "cell_extent",
        "cell_feret", "cell_major_minor", "cell_bbox_hull",
        "cell_moments", "cell_equiv_diam",
    ]),
    ("group_shape_2d_soma", [
        "soma_area", "soma_perimeter", "soma_solidity", "soma_elongation",
        "soma_form_factor", "soma_roundness", "soma_eccentricity",
        "soma_convexity",
    ]),
    ("group_shape_2d_nucleus", [
        "nucleus_area", "nucleus_perimeter", "nucleus_solidity",
        "nucleus_elongation", "nucleus_form_factor", "nucleus_roundness",
        "nucleus_eccentricity", "nucleus_convexity",
    ]),
    ("group_shape_3d_cell", [
        "cell_volume", "cell_surface", "cell_sphericity", "cell_sav",
        "cell_convexity_3d", "cell_extent_3d", "cell_axes_3d",
    ]),
    ("group_shape_3d_soma", [
        "soma_volume", "soma_elongation_3d", "soma_flatness",
        "soma_extent_3d", "soma_convexity_3d",
    ]),
    ("group_shape_3d_nucleus", [
        "nucleus_volume", "nucleus_elongation_3d", "nucleus_flatness",
        "nucleus_extent_3d", "nucleus_convexity_3d",
    ]),
    ("group_intensity_cell", [
        "cell_mean_intensity", "cell_intensity_stats",
        "cell_intensity_percentiles", "cell_intensity_entropy_gini",
        "cell_intensity_highlow", "cell_snr_cnr_sbr",
    ]),
    ("group_intensity_soma", [
        "soma_mean_intensity", "soma_intensity_stats",
        "soma_intensity_entropy",
    ]),
    ("group_intensity_nucleus", [
        "nucleus_mean_intensity", "nucleus_intensity_stats",
        "nucleus_intensity_entropy",
    ]),
    ("group_texture_cell", [
        "cell_local_variance", "cell_glcm", "cell_fractal",
    ]),
    ("group_texture_soma", [
        "soma_glcm",
    ]),
    ("group_texture_nucleus", [
        "nucleus_glcm",
    ]),
    ("group_nucleus", [
        "cell_nuc_cyto", "cell_nuc_intensity", "cell_centroid_dist",
    ]),
    ("group_skeleton", [
        "cell_process", "cell_skeleton",
    ]),
    ("group_granules", [
        "granule_intensity", "granule_size_count", "granule_shape",
        "granule_spatial", "granule_sbr",
    ]),
]

def default_metrics_config(is_2d=None):
    cfg = {}
    for _grp, items in METRIC_GROUPS:
        for it in items:
            cfg[it] = False

    basic_2d = [
        "cell_area", "cell_perimeter", "cell_solidity",
        "cell_elongation", "cell_form_factor", "cell_roundness",
        "cell_aspect_ratio", "cell_eccentricity",
        "soma_area", "soma_perimeter", "soma_solidity",
        "soma_elongation", "soma_form_factor", "soma_roundness",
        "soma_eccentricity",
        "nucleus_area", "nucleus_perimeter", "nucleus_solidity",
        "nucleus_elongation", "nucleus_form_factor", "nucleus_roundness",
        "nucleus_eccentricity",
    ]
    basic_3d = [
        "cell_volume", "cell_surface", "cell_sphericity",
        "soma_volume", "nucleus_volume",
        "cell_nuc_cyto",
    ]

    if is_2d is None:
        basic_metrics = basic_2d + basic_3d
    elif is_2d:
        basic_metrics = basic_2d + ["cell_nuc_cyto"]
    else:
        basic_metrics = basic_3d

    for it in basic_metrics:
        if it in cfg:
            cfg[it] = True

    return cfg


def metric_enabled(cfg, name, is_2d=None):
    if not cfg:
        return True
    if is_2d is not None:
        if is_2d and name in _3D_ONLY_METRICS:
            return False
        if not is_2d and name in _2D_ONLY_METRICS:
            return False
    return bool(cfg.get(name, True))

class AnalysisContext:
    def __init__(self):
        self.astrocyte_channel = 0
        self.dapi_channel = 1
        self.cell_type = "astrocytes"
        self.is_2d = False
        self.bridge = None
        self.app = None
        self.lang = "ru"
        self.image_mode = "confocal"
        self._soma_labels = None
        self._channels_initialized = False
        self._optics_initialized = False
        self.downsample = 1
        self.min_cell_volume = 100
        self.xy_scale_um = None
        self.z_scale_um = None
        self.scale_initialized = False
        self.deconvolve_images = True
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
        self._granule_deconv_params = {
            'na': 0.95,
            'emission_wavelength_um': 0.520,
            'wiener_epsilon': 0.005
        }
        self._current_base_name = None
        self._current_out_dir = None
        self._growth_mode = "normal"
        self._gui_scale_mode = "auto"
        self._gui_deconv_mode = "auto"
        self._gui_do_filter = True
        self._gui_settings_locked = False
        self.has_dapi = True
        self.has_granules = False
        self.auto_soma = False
        self.n_channels_original = 1
        self.cell_channel_idx = None
        self.dapi_channel_idx = None
        self.granule_channel_idx = None
        self.selected_channels = None
        self._preprocess_decided = False
        self.crop_box = None
        self.allow_downsample = True
        self.metrics_config = default_metrics_config()

    def t(self, key, **kwargs):
        if getattr(self, "app", None) is not None:
            lang = self.app.var_lang.get()
        else:
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

            selected_idx = [cell_idx]

            if self.ctx.has_dapi:
                dapi_idx = self.ctx.dapi_channel_idx
                if (dapi_idx is None or dapi_idx < 0 or
                        dapi_idx >= n_ch or dapi_idx == cell_idx):
                    alt = [k for k in range(n_ch) if k != cell_idx and k not in selected_idx]
                    if alt:
                        dapi_idx = alt[0]
                    else:
                        self.ctx.has_dapi = False
                        dapi_idx = None
                if self.ctx.has_dapi and dapi_idx not in selected_idx:
                    selected_idx.append(dapi_idx)

            if getattr(self.ctx, "has_granules", False):
                gran_idx = getattr(self.ctx, "granule_channel_idx", None)
                if (gran_idx is None or gran_idx < 0 or
                        gran_idx >= n_ch or gran_idx in selected_idx):
                    alt = [k for k in range(n_ch) if k not in selected_idx]
                    if alt:
                        gran_idx = alt[0]
                    else:
                        self.ctx.has_granules = False
                        gran_idx = None
                if self.ctx.has_granules and gran_idx not in selected_idx:
                    selected_idx.append(gran_idx)

            self.ctx.selected_channels = selected_idx

            img = img[selected_idx, ...]

            if self.ctx.has_dapi and img.shape[0] < 2:
                self.ctx.log(self.ctx.t("log.warn_only_one_channel", n=img.shape[0]))
                img = np.concatenate([img, img], axis=0)

            c, z, y, x = img.shape

            crop_box = getattr(self.ctx, "crop_box", None)
            if crop_box is not None:
                cy0, cy1, cx0, cx1 = crop_box
                cy0 = max(0, min(int(cy0), y))
                cy1 = max(cy0 + 1, min(int(cy1), y))
                cx0 = max(0, min(int(cx0), x))
                cx1 = max(cx0 + 1, min(int(cx1), x))
                img = img[:, :, cy0:cy1, cx0:cx1]
                c, z, y, x = img.shape

            self.ctx.is_2d = (z == 1)
            if self.ctx.is_2d:
                self.ctx.log(self.ctx.t("log.2d_image"))
            else:
                self.ctx.log(self.ctx.t("log.3d_stack", n=z))

            if not getattr(self.ctx, "allow_downsample", True):
                return img.astype(np.float32)

            ds = max(1, int(self.ctx.downsample))
            if ds == 1:
                return img.astype(np.float32)

            ny, nx = y // ds, x // ds
            if ny < 1 or nx < 1:
                return img.astype(np.float32)

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
            return

        if meta['source'] == 'none' or (meta['xy_um'] is None and meta['z_um'] is None):
            self.ctx.log(self.ctx.t("log.no_voxel_metadata"))
            self.ctx.log(self.ctx.t("log.use_panel_settings"))
            if self.ctx.xy_scale_um is None:
                self.ctx.xy_scale_um = 0.144
            if self.ctx.z_scale_um is None:
                self.ctx.z_scale_um = 1.0
            self.ctx.scale_initialized = True
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
        self.ctx.log(self.ctx.t("log.scale_value",
                                xy=self.ctx.xy_scale_um,
                                z=self.ctx.z_scale_um))

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

        supported_exts = ['.nd2', '.oir', '.czi', '.oib', '.oif', '.lif', '.tif', '.tiff']
        if ext not in supported_exts:
            return result

        try:
            from bioio import BioImage
            import bioio_openreadout

            img = BioImage(filepath, reader=bioio_openreadout.Reader)

            na_global = None
            channels = []

            if hasattr(img, 'ome_metadata') and img.ome_metadata is not None:
                ome = img.ome_metadata

                for instrument in getattr(ome, 'instruments', []):
                    for objective in getattr(instrument, 'objectives', []):
                        if objective.lens_na is not None:
                            na_global = float(objective.lens_na)
                            break
                    if na_global is not None:
                        break

                for image in getattr(ome, 'images', []):
                    pixels = getattr(image, 'pixels', None)
                    if pixels is None:
                        continue
                    for ch in getattr(pixels, 'channels', []):
                        em_um = None
                        if ch.emission_wavelength is not None:
                            em_um = float(ch.emission_wavelength) / 1000.0
                        channels.append({
                            'index': len(channels),
                            'name': ch.name or ch.id or f"Channel {len(channels)}",
                            'emission_um': em_um,
                            'na': na_global,
                            'magnification': None
                        })
                    break

            if not channels:
                names = getattr(img, 'channel_names', [])
                for i, name in enumerate(names):
                    channels.append({
                        'index': i,
                        'name': name,
                        'emission_um': None,
                        'na': na_global,
                        'magnification': None
                    })

            if not channels:
                for i in range(2):
                    channels.append({
                        'index': i,
                        'name': f"Channel {i}",
                        'emission_um': None,
                        'na': na_global,
                        'magnification': None
                    })

            result['channels'] = channels
            self._assign_dapi_by_metadata(result)
            result['source'] = 'bioio_openreadout'

        except ImportError:
            result['source'] = 'bioio_not_installed'
        except Exception as e:
            result['source'] = f'bioio_error: {str(e)}'

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
                                            ex=0.0))
                else:
                    self.ctx.log(self.ctx.t("log.channel_line_no_name",
                                            idx=ch['index'],
                                            em=(ch.get('emission_um') or 0) * 1000,
                                            ex=0.0))

        auto_cell_idx = optics.get('astro_channel_index')
        auto_dapi_idx = optics.get('dapi_channel_index')

        img = self.file_io._read_image_universal(filepath)
        if img is None:
            self.ctx.log(self.ctx.t("log.file_read_failed"))
            self.ctx.has_dapi = False
            self.ctx.has_granules = False
            self.ctx.cell_channel_idx = 0
            self.ctx.dapi_channel_idx = None
            self.ctx.granule_channel_idx = None
            self.ctx.selected_channels = [0]
            self.ctx.astrocyte_channel = 0
            self.ctx.dapi_channel = 0
            self.ctx.metrics_config = default_metrics_config()
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
            self.ctx.has_granules = False
            self.ctx.cell_channel_idx = 0
            self.ctx.dapi_channel_idx = None
            self.ctx.granule_channel_idx = None
            self.ctx.selected_channels = [0]
            self.ctx.astrocyte_channel = 0
            self.ctx.dapi_channel = 0
            self.ctx.metrics_config = default_metrics_config()
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
            self.ctx.has_granules = False
            self.ctx.granule_channel_idx = None
            if has_dapi_default:
                self.ctx.selected_channels = [cell_default, dapi_default]
            else:
                self.ctx.selected_channels = [cell_default]
            self.ctx.metrics_config = default_metrics_config()
        else:
            self.ctx.cell_channel_idx = int(res.get("cell_idx", cell_default))
            self.ctx.has_dapi = bool(res.get("has_dapi", has_dapi_default))

            d_idx = res.get("dapi_idx")
            self.ctx.dapi_channel_idx = int(d_idx) if d_idx is not None else None

            self.ctx.has_granules = bool(res.get("has_granules", False))

            g_idx = res.get("granule_idx")
            self.ctx.granule_channel_idx = int(g_idx) if g_idx is not None else None

            sel = [self.ctx.cell_channel_idx]
            if self.ctx.has_dapi and self.ctx.dapi_channel_idx is not None:
                sel.append(self.ctx.dapi_channel_idx)
            if self.ctx.has_granules and self.ctx.granule_channel_idx is not None:
                sel.append(self.ctx.granule_channel_idx)
            self.ctx.selected_channels = sorted(list(set(sel)))

            mc = res.get("metrics_config", None)
            if isinstance(mc, dict) and mc:
                self.ctx.metrics_config = dict(mc)
            else:
                self.ctx.metrics_config = default_metrics_config()

        self.ctx.astrocyte_channel = self.ctx.cell_channel_idx
        self.ctx.dapi_channel = self.ctx.dapi_channel_idx if self.ctx.has_dapi else 0

        self.ctx.log(self.ctx.t("log.cells_channel", idx=self.ctx.cell_channel_idx))
        if self.ctx.has_dapi:
            self.ctx.log(self.ctx.t("log.dapi_channel", idx=self.ctx.dapi_channel_idx))
        else:
            self.ctx.log(self.ctx.t("log.dapi_not_used"))
        if self.ctx.has_granules:
            self.ctx.log(self.ctx.t("log.granule_channel", idx=self.ctx.granule_channel_idx))
        self.ctx.log(self.ctx.t("log.channels_in_analysis",
                                ch=self.ctx.selected_channels))

    def extract_and_confirm_optics(self, filepath):
        self.ctx.log(self.ctx.t("log.optics_title"))
        self.ctx.log(self.ctx.t("log.file_label", name=os.path.basename(filepath)))

        optics = self._extract_optics_from_file(filepath)
        channels = optics.get('channels', [])

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

        dapi_idx = self.ctx.dapi_channel_idx
        astro_idx = self.ctx.cell_channel_idx

        astro_ch = channels[astro_idx] if (astro_idx is not None and astro_idx < len(channels)) else {}
        dapi_ch = channels[dapi_idx] if (dapi_idx is not None and dapi_idx < len(channels)) else {}

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

    def _deconvolve_3d_stack(self, image_stack, channel_label_key="channel.cells"):
        self.ctx.log(self.ctx.t("log.deconv_start_channel", ch=self.ctx.t(channel_label_key)))

        z_dim, y_dim, x_dim = image_stack.shape

        sigma_xy_bg = 2.0
        sigma_z_bg = 5.0
        radius_xy_bg = int(sigma_xy_bg * 3) + 1
        radius_z_bg = int(sigma_z_bg * 3) + 1

        z = np.arange(-radius_z_bg, radius_z_bg + 1)
        y = np.arange(-radius_xy_bg, radius_xy_bg + 1)
        x = np.arange(-radius_xy_bg, radius_xy_bg + 1)
        Z, Y, X = np.meshgrid(z, y, x, indexing='ij')

        psf_kernel_bg = np.exp(-((X**2 + Y**2) / (2 * sigma_xy_bg**2) +
                                 (Z**2) / (2 * sigma_z_bg**2)))
        psf_kernel_bg /= np.sum(psf_kernel_bg)

        img_bg = convolve(image_stack, psf_kernel_bg, mode='reflect')
        img_subtracted = image_stack - img_bg * 1.0
        img_subtracted = np.clip(img_subtracted, 0, None)

        sigma_unsharp = 0.3
        img_blur_small = gaussian_filter(img_subtracted, sigma=sigma_unsharp)
        unsharp_weight = 0.4
        img_sharp = img_subtracted + unsharp_weight * (img_subtracted - img_blur_small)
        img_sharp = np.clip(img_sharp, 0, None)

        img_deconv = median_filter(img_sharp, size=(1, 2, 2), mode='reflect')

        original_sum = np.sum(image_stack)
        if original_sum > 0 and np.sum(img_deconv) > 0:
            scale_factor = original_sum / (np.sum(img_deconv) + 1e-6)
            img_deconv = img_deconv * scale_factor

        to_delete = self.ctx.bridge._ask(
            "slice_selector",
            stack=img_deconv,
            title=self.ctx.t("dlg.slices.title_selector_channel", ch=self.ctx.t(channel_label_key)),
        )

        if to_delete and len(to_delete) > 0:
            self.ctx.log(self.ctx.t("log.deconv_slices_deleted_channel",
                                    ch=self.ctx.t(channel_label_key), n=len(to_delete), z=to_delete))
            img_deconv[to_delete, :, :] = 0.0
        else:
            self.ctx.log(self.ctx.t("log.deconv_all_slices_kept_channel", ch=self.ctx.t(channel_label_key)))

        self.ctx.log(self.ctx.t("log.deconv_done_channel", ch=self.ctx.t(channel_label_key), shape=img_deconv.shape))
        return img_deconv.astype(np.float32)

    def _deconvolve_3d_stack_dapi(self, image_stack, pixel_xy_um, pixel_z_um,
                                  channel_label_key="channel.dapi"):
        self.ctx.log(self.ctx.t("log.deconv_start_channel", ch=self.ctx.t(channel_label_key)))

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
            title=self.ctx.t("dlg.slices.title_selector_channel", ch=self.ctx.t(channel_label_key)),
        )

        if to_delete and len(to_delete) > 0:
            self.ctx.log(self.ctx.t("log.deconv_slices_deleted_channel",
                                    ch=self.ctx.t(channel_label_key), n=len(to_delete), z=to_delete))
            img_deconv[to_delete, :, :] = 0.0
        else:
            self.ctx.log(self.ctx.t("log.deconv_all_slices_kept_channel", ch=self.ctx.t(channel_label_key)))

        self.ctx.log(self.ctx.t("log.deconv_done_channel", ch=self.ctx.t(channel_label_key), shape=img_deconv.shape))
        return img_deconv.astype(np.float32)

class Segmentation:
    def __init__(self, ctx):
        self.ctx = ctx

    def _compute_astro_norm(self, astrocyte_stack, image_mode="confocal", is_granules=False):
        if is_granules:
            self.ctx.log(self.ctx.t("log.granules_title"))

        z_dim, y_dim, x_dim = astrocyte_stack.shape
        astro = astrocyte_stack.astype(np.float32)

        if image_mode == "epifluorescence":
            astro_norm = self._preprocess_epi_log(astro)
            return astro_norm

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
        mean_structure = np.zeros_like(astro)
        std_structure = np.zeros_like(astro)
        for smap in structure_maps:
            max_structure = np.maximum(max_structure, smap)
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

        mip = np.max(astro, axis=0)
        sigma_bg = max(30, min(y_dim, x_dim) // 15)
        mip_bg = gaussian_filter(mip, sigma=sigma_bg)
        mip_bg_safe = np.maximum(mip_bg, 1e-6)
        local_contrast = mip / mip_bg_safe
        contrast_conf = np.clip((local_contrast - 1.0) / 1.0, 0.0, 1.0).astype(np.float32)

        nz_mip = mip[mip > 0]
        if nz_mip.size > 0:
            p50, p95 = np.percentile(nz_mip, [50, 95])
            bright_conf = np.clip((mip - p50) / (p95 - p50 + 1e-6), 0.0, 1.0).astype(np.float32)
        else:
            bright_conf = np.zeros_like(mip, dtype=np.float32)

        signal_conf = np.clip(0.6 * contrast_conf + 0.4 * bright_conf, 0.0, 1.0)
        signal_conf = gaussian_filter(signal_conf, sigma=2.0)

        grain_mip = np.max(grain_index, axis=0)
        nz_g = grain_mip[grain_mip > 0]
        if nz_g.size > 0:
            g_med = float(np.median(nz_g))
            g_p90 = float(np.percentile(nz_g, 90))
        else:
            g_med = 0.0
            g_p90 = 0.0
        g_thr = max(g_med, g_p90 * 0.6)
        edge_mask_mip = (grain_mip >= g_thr) & (signal_conf > 0.35)
        edge_mask_mip = edge_mask_mip.astype(np.float32)

        radius = max(6, min(y_dim, x_dim) // 12)
        ksize = 2 * radius + 1
        yy, xx = np.mgrid[-radius:radius + 1, -radius:radius + 1]
        dist2 = yy * yy + xx * xx
        ring_lo = (radius * 0.5) ** 2
        ring_hi = float(radius * radius)
        ring_kernel = ((dist2 >= ring_lo) & (dist2 <= ring_hi)).astype(np.float32)
        ring_kernel /= (ring_kernel.sum() + 1e-6)

        directional_hits = np.zeros((8, y_dim, x_dim), dtype=np.float32)
        directions = [
            (0, 1), (1, 1), (1, 0), (1, -1),
            (0, -1), (-1, -1), (-1, 0), (-1, 1),
        ]
        probe_len = radius
        for di, (dy, dx) in enumerate(directions):
            hit = np.zeros((y_dim, x_dim), dtype=np.float32)
            for step in range(1, probe_len + 1):
                sy = dy * step
                sx = dx * step
                shifted = np.zeros_like(edge_mask_mip)
                ys0 = max(0, -sy); ys1 = y_dim - max(0, sy)
                xs0 = max(0, -sx); xs1 = x_dim - max(0, sx)
                if ys1 <= ys0 or xs1 <= xs0:
                    continue
                shifted[ys0 + sy:ys1 + sy, xs0 + sx:xs1 + sx] = edge_mask_mip[ys0:ys1, xs0:xs1]
                hit = np.maximum(hit, shifted)
                if step >= 2 and np.all(hit > 0):
                    break
            directional_hits[di] = hit

        surround_count = directional_hits.sum(axis=0)
        surround_conf = np.clip(surround_count / 4.0, 0.0, 1.0).astype(np.float32)
        surround_conf = gaussian_filter(surround_conf, sigma=1.5)

        weight_map_mip = signal_conf * (0.25 + 0.75 * surround_conf)
        weight_map_mip = np.clip(weight_map_mip, 0.0, 1.0).astype(np.float32)
        weight_map_mip = gaussian_filter(weight_map_mip, sigma=(1.0, 1.0))

        astro_clean = np.zeros_like(astro, dtype=np.float32)
        for z in range(z_dim):
            slice_data = astro[z]
            slice_weight = weight_map_mip

            base_clean = slice_data - mip_bg * 0.7 * (1.0 - slice_weight)

            if np.max(base_clean) > 0:
                slice_bg = gaussian_filter(
                    base_clean,
                    sigma=max(20, min(y_dim, x_dim) // 20)
                )
                corrected = base_clean - slice_bg * 0.7 * (1.0 - slice_weight)
            else:
                corrected = base_clean

            astro_clean[z] = np.clip(corrected, 0, None)

        astro_clean = median_filter(astro_clean, size=(1, 2, 2))

        if np.max(astro_clean) == 0:
            return None

        p2, p98 = np.percentile(astro_clean[astro_clean > 0], [2, 98])
        astro_norm = np.clip((astro_clean - p2) / (p98 - p2 + 1e-6), 0, 1)

        self.ctx._weight_map = weight_map_mip[np.newaxis, :, :].repeat(z_dim, axis=0).astype(np.float32)
        self.ctx._grain_index_map = grain_index.astype(np.float32)

        return astro_norm

    def segment_nuclei_local_from_soma(self, soma_mask, dapi_stack):
        self.ctx.log(self.ctx.t("log.nuclei_title"))

        if not getattr(self.ctx, "has_dapi", True):
            self.ctx.log(self.ctx.t("log.nuclei_no_dapi"))
            return np.zeros_like(dapi_stack, dtype=np.int32)

        from scipy.ndimage import (gaussian_filter, binary_fill_holes,
                                   label as nd_label, binary_dilation)
        from skimage.morphology import remove_small_objects
        from skimage.filters import threshold_triangle
        from skimage.measure import find_contours

        z_dim, y_dim, x_dim = dapi_stack.shape
        clean_nuclei = np.zeros_like(dapi_stack, dtype=np.int32)
        unique_ids = np.unique(soma_mask)[1:]
        if len(unique_ids) == 0:
            self.ctx.log(self.ctx.t("log.nuclei_no_soma"))
            return clean_nuclei

        xy_um = float(self.ctx.xy_scale_um or 1.0) * float(self.ctx.downsample)
        z_um = float(self.ctx.z_scale_um or 1.0)

        dapi_smooth = gaussian_filter(
            dapi_stack.astype(np.float32), sigma=(0.8, 1.0, 1.0)
        )

        nz = dapi_smooth[dapi_smooth > 0]
        if nz.size < 100:
            return clean_nuclei

        try:
            signal_thr = float(threshold_triangle(nz))
        except Exception:
            signal_thr = float(np.percentile(nz, 35))
        p20 = float(np.percentile(nz, 20))
        p99 = float(np.percentile(nz, 99))
        dyn_range = max(p99 - p20, 1e-6)
        signal_thr = min(signal_thr, p20 + 0.6 * dyn_range)
        signal_thr = max(signal_thr, p20 + 0.15 * dyn_range)

        n_steps = 25
        thresholds = np.linspace(p99, signal_thr, n_steps)
        persistence_map = np.zeros_like(dapi_smooth, dtype=np.float32)

        prev_labels = np.zeros_like(dapi_smooth, dtype=np.int32)
        prev_n = 0
        birth_thr = {}

        for step, thr in enumerate(thresholds):
            mask = dapi_smooth >= thr
            mask = remove_small_objects(mask, min_size=10)
            lbl, n = nd_label(mask, structure=np.ones((3, 3, 3)))

            if n == 0:
                prev_labels = np.zeros_like(dapi_smooth, dtype=np.int32)
                prev_n = 0
                continue

            if prev_n == 0:
                for k in range(1, n + 1):
                    birth_thr[k] = thr
                prev_labels = lbl
                prev_n = n
                continue

            for k in range(1, n + 1):
                comp_mask = (lbl == k)
                overlap = prev_labels[comp_mask]
                overlap = overlap[overlap > 0]
                if len(overlap) == 0:
                    birth_thr[k] = thr
                    continue
                vals, counts = np.unique(overlap, return_counts=True)
                best_prev = int(vals[np.argmax(counts)])
                birth_thr[k] = birth_thr.get(best_prev, thr)

            prev_labels = lbl
            prev_n = n

        for step, thr in enumerate(thresholds):
            mask = dapi_smooth >= thr
            mask = remove_small_objects(mask, min_size=10)
            lbl, n = nd_label(mask, structure=np.ones((3, 3, 3)))
            if n == 0:
                continue
            for k in range(1, n + 1):
                if k not in birth_thr:
                    continue
                lifespan = birth_thr[k] - thr
                if lifespan > 0:
                    comp_vals = lbl == k
                    current_max = persistence_map[comp_vals].max()
                    if lifespan > current_max:
                        persistence_map[comp_vals] = lifespan

        persistence_thr = 0.15 * dyn_range
        persistent_mask = persistence_map >= persistence_thr
        persistent_mask = remove_small_objects(persistent_mask, min_size=10)
        persistent_mask = binary_fill_holes(persistent_mask)

        lbl_p, n_p = nd_label(persistent_mask, structure=np.ones((3, 3, 3)))

        seed_info = []
        for cid in unique_ids:
            soma_cid = (soma_mask == cid)
            if not np.any(soma_cid):
                continue
            coords = np.argwhere(soma_cid)
            center = coords.mean(axis=0)
            cz = int(round(center[0]))
            cy = int(round(center[1]))
            cx = int(round(center[2]))
            cz = max(0, min(z_dim - 1, cz))
            cy = max(0, min(y_dim - 1, cy))
            cx = max(0, min(x_dim - 1, cx))

            if not persistent_mask[cz, cy, cx]:
                search_r = 4
                z_lo, z_hi = max(0, cz - search_r), min(z_dim, cz + search_r + 1)
                y_lo, y_hi = max(0, cy - search_r), min(y_dim, cy + search_r + 1)
                x_lo, x_hi = max(0, cx - search_r), min(x_dim, cx + search_r + 1)
                sub = dapi_smooth[z_lo:z_hi, y_lo:y_hi, x_lo:x_hi]
                sub_p = persistent_mask[z_lo:z_hi, y_lo:y_hi, x_lo:x_hi]
                if np.any(sub_p):
                    masked = np.where(sub_p, sub, -np.inf)
                    idx = np.unravel_index(np.argmax(masked), masked.shape)
                    cz = int(idx[0] + z_lo)
                    cy = int(idx[1] + y_lo)
                    cx = int(idx[2] + x_lo)
                else:
                    continue

            comp_id = int(lbl_p[cz, cy, cx])
            if comp_id == 0:
                continue

            seed_info.append({
                'cid': int(cid),
                'seed': (cz, cy, cx),
                'comp_id': comp_id
            })

        if not seed_info:
            self.ctx.log(self.ctx.t("log.nuclei_no_soma"))
            return clean_nuclei

        comp_to_seeds = {}
        for info in seed_info:
            comp_to_seeds.setdefault(info['comp_id'], []).append(info)

        for comp_id, infos in comp_to_seeds.items():
            comp_mask = (lbl_p == comp_id)

            cid = infos[0]['cid']
            soma_here = (soma_mask == cid)
            soma_dilated = binary_dilation(soma_here, structure=np.ones((3, 3, 3)), iterations=1)

            if len(infos) == 1:
                nuc_mask = comp_mask & soma_dilated
                nuc_mask = binary_fill_holes(nuc_mask)
                nuc_mask = remove_small_objects(nuc_mask, min_size=10)
                if np.sum(nuc_mask) >= 10:
                    clean_nuclei[nuc_mask] = cid
                continue

            cz, cy, cx = infos[0]['seed']
            slice_mask = (comp_mask & soma_dilated)[cz]
            if not np.any(slice_mask):
                for info in infos:
                    c_info = info['cid']
                    c_soma = binary_dilation(soma_mask == c_info, structure=np.ones((3, 3, 3)), iterations=1)
                    c_mask = comp_mask & c_soma
                    if np.sum(c_mask) >= 10:
                        clean_nuclei[c_mask] = c_info
                continue

            contours = find_contours(slice_mask.astype(float), 0.5)
            if not contours:
                for info in infos:
                    c_info = info['cid']
                    c_soma = binary_dilation(soma_mask == c_info, structure=np.ones((3, 3, 3)), iterations=1)
                    c_mask = comp_mask & c_soma
                    if np.sum(c_mask) >= 10:
                        clean_nuclei[c_mask] = c_info
                continue

            contour = max(contours, key=len)
            n_pts = len(contour)
            if n_pts < 20:
                for info in infos:
                    c_info = info['cid']
                    c_soma = binary_dilation(soma_mask == c_info, structure=np.ones((3, 3, 3)), iterations=1)
                    c_mask = comp_mask & c_soma
                    if np.sum(c_mask) >= 10:
                        clean_nuclei[c_mask] = c_info
                continue

            concavity_scores = np.zeros(n_pts, dtype=np.float32)
            window = max(3, n_pts // 20)
            for i in range(n_pts):
                p_prev = contour[(i - window) % n_pts]
                p_curr = contour[i]
                p_next = contour[(i + window) % n_pts]
                v1 = p_curr - p_prev
                v2 = p_next - p_curr
                cross = v1[0] * v2[1] - v1[1] * v2[0]
                concavity_scores[i] = cross

            candidate_idx = []
            for i in range(n_pts):
                if concavity_scores[i] < np.percentile(concavity_scores, 15):
                    candidate_idx.append(i)

            split_points = []
            min_dist = max(3, n_pts // 10)
            for idx in candidate_idx:
                if not split_points:
                    split_points.append(idx)
                else:
                    if idx - split_points[-1] > min_dist and (n_pts - idx + split_points[0]) > min_dist:
                        split_points.append(idx)

            if len(split_points) < 2:
                for info in infos:
                    c_info = info['cid']
                    c_soma = binary_dilation(soma_mask == c_info, structure=np.ones((3, 3, 3)), iterations=1)
                    c_mask = comp_mask & c_soma
                    if np.sum(c_mask) >= 10:
                        clean_nuclei[c_mask] = c_info
                continue

            seed_slice_coords = []
            for info in infos:
                sz, sy, sx = info['seed']
                seed_slice_coords.append((info['cid'], sy, sx))

            assigned = np.zeros_like(slice_mask, dtype=np.int32)
            for c_info, sy, sx in seed_slice_coords:
                if 0 <= sy < slice_mask.shape[0] and 0 <= sx < slice_mask.shape[1]:
                    seed_pt = np.zeros_like(slice_mask, dtype=bool)
                    seed_pt[sy, sx] = True
                    seed_pt = binary_dilation(seed_pt, structure=np.ones((3, 3)), iterations=2)
                    seed_pt = seed_pt & slice_mask
                    assigned[seed_pt] = c_info

            queue = deque()
            for c_info, sy, sx in seed_slice_coords:
                if assigned[sy, sx] == c_info:
                    queue.append((sy, sx))

            while queue:
                y, x = queue.popleft()
                c_info = assigned[y, x]
                for dy, dx in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1)]:
                    ny, nx = y + dy, x + dx
                    if 0 <= ny < slice_mask.shape[0] and 0 <= nx < slice_mask.shape[1]:
                        if slice_mask[ny, nx] and assigned[ny, nx] == 0:
                            assigned[ny, nx] = c_info
                            queue.append((ny, nx))

            for info in infos:
                c_info = info['cid']
                z_slice = info['seed'][0]
                slice_assigned = assigned == c_info
                if not np.any(slice_assigned):
                    continue

                c_soma = binary_dilation(soma_mask == c_info, structure=np.ones((3, 3, 3)), iterations=1)
                c_mask = comp_mask & c_soma

                current_nuc = c_mask.copy()
                current_nuc[z_slice] = slice_assigned & slice_mask
                current_nuc = binary_fill_holes(current_nuc)
                current_nuc = remove_small_objects(current_nuc, min_size=10)

                if np.sum(current_nuc) >= 10:
                    clean_nuclei[current_nuc] = c_info

        self.ctx.log(self.ctx.t("log.nuclei_found",
                                n=len(np.unique(clean_nuclei)) - 1))
        return clean_nuclei

    def _grow_cells_from_soma(self, astro_norm, soma_mask):
        self.ctx.log(self.ctx.t("log.cells_title"))

        if astro_norm is None:
            self.ctx.log(self.ctx.t("log.cells_no_signal"))
            return np.zeros_like(soma_mask, dtype=np.int32)

        unique_ids = np.unique(soma_mask)[1:]
        if len(unique_ids) == 0:
            self.ctx.log(self.ctx.t("log.cells_no_soma"))
            return np.zeros_like(soma_mask, dtype=np.int32)

        return self._grow_cells_from_somas(astro_norm, soma_mask, unique_ids).astype(np.int32)

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

    def _grow_cells_by_components(self, astro_norm, soma_mask, growth_mask,
                                  unique_ids, dip_hard_threshold=0.5):
        import heapq
        from scipy.ndimage import laplace

        z_dim, y_dim, x_dim = astro_norm.shape
        astro = astro_norm.astype(np.float32)

        is_soft = getattr(self.ctx, "_growth_mode", "normal") == "soft"

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

        if is_soft:
            dip_hard = dip_norm > 0.95
        else:
            dip_hard = dip_norm > dip_hard_threshold

        weight[dip_hard] = 0.0
        weight = np.clip(weight, 0.0, None).astype(np.float32)

        weight_in_graph = weight[graph_mask]
        if len(weight_in_graph) > 0:
            weight_median = float(np.median(weight_in_graph))
        else:
            weight_median = 0.0

        if is_soft:
            forbidden_mask = np.zeros_like(soma_bool)
        else:
            forbidden_mask = weight < (weight_median * 0.3)
            forbidden_mask = forbidden_mask & graph_mask
            forbidden_mask = forbidden_mask & (~soma_bool)

        eps = 1e-6
        z_step_penalty = 2.0

        straight_penalty_factor = 3.0

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

            dir_active = np.full((n_active, 3), np.nan, dtype=np.float32)

            heap = []
            for i, (zz, yy, xx) in enumerate(coords_list):
                if soma_bool[zz, yy, xx]:
                    dist_active[i] = 0.0
                    source_active[i] = int(soma_mask[zz, yy, xx])
                    dir_active[i] = [0.0, 0.0, 0.0]
                    heapq.heappush(heap, (0.0, i))

            while heap:
                d, i = heapq.heappop(heap)
                if d > dist_active[i]:
                    continue

                z, y, x = coords_list[i]
                w_here = weight[z, y, x]
                src_here = source_active[i]
                dir_here = dir_active[i]

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

                    step_vec = np.array([dz, dy, dx], dtype=np.float32)
                    step_vec /= np.linalg.norm(step_vec)

                    if not np.isnan(dir_here).any():
                        dot = np.dot(dir_here, step_vec)
                        if dot < 0.5:
                            penalty = 1.0 + (1.0 - dot) * straight_penalty_factor
                            step_cost *= penalty

                    nd = d + step_cost

                    if nd < dist_active[j]:
                        dist_active[j] = nd
                        source_active[j] = src_here
                        dir_active[j] = step_vec
                        heapq.heappush(heap, (nd, j))

            for i, (zz, yy, xx) in enumerate(coords_list):
                src_id = int(source_active[i])
                if src_id > 0:
                    final_cells[zz, yy, xx] = src_id

        return final_cells

    def _preprocess_epi_log(self, astro):
        from scipy.ndimage import laplace

        sigma_xy = 1.0
        sigma_z = 1.0
        blurred = gaussian_filter(astro, sigma=(sigma_z, sigma_xy, sigma_xy), mode='reflect')
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

        lap_coh = gaussian_filter(log_map, sigma=(0.0, sigma_s, sigma_s), mode='reflect')
        base_level = gaussian_filter(np.abs(lap_coh), sigma=(0.0, sigma_l, sigma_l), mode='reflect')
        mask = np.abs(lap_coh) > (k * base_level)

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

    def _multilevel_soma_threshold(self, mip):
        img = mip.astype(np.float32)
        nz = img[img > 0]
        if nz.size < 50:
            return None, None, None
        try:
            classes = 5
            thresholds = threshold_multiotsu(nz, classes=classes)
        except Exception:
            return None, None, None
        if thresholds is None or len(thresholds) == 0:
            return None, None, None
        thresholds = np.sort(np.asarray(thresholds, dtype=np.float32))

        min_area_px = 30
        n_list = []
        for t in thresholds:
            bin_mask = (img >= float(t))
            labeled, n_obj = label(bin_mask, structure=np.ones((3, 3)))
            if n_obj == 0:
                n_list.append(0)
                continue
            sizes = np.bincount(labeled.ravel())
            if sizes.size > 1:
                sizes = sizes[1:]
                n_ok = int(np.sum(sizes >= min_area_px))
            else:
                n_ok = 0
            n_list.append(n_ok)

        n_arr = np.asarray(n_list, dtype=np.float32)

        growing_idx = []
        for i in range(len(n_arr) - 1):
            denom = n_arr[i + 1] + n_arr[i]
            if denom <= 0:
                continue
            c = (n_arr[i + 1] - n_arr[i]) / denom
            if c >= 0:
                growing_idx.append((i, c))

        if not growing_idx:
            t_s_idx = len(thresholds) // 2
        else:
            c_vals = np.array([c for _, c in growing_idx], dtype=np.float32)
            c_min = float(np.min(c_vals))
            cand = [growing_idx[j][0] for j in range(len(growing_idx))
                    if abs(c_vals[j] - c_min) < 1e-6]

            if len(cand) == 1:
                t_s_idx = cand[0]
            else:
                best_idx = cand[0]
                best_delta = float("inf")
                for idx in cand:
                    if idx + 1 < len(c_vals):
                        delta = abs(float(c_vals[idx + 1] - c_vals[idx]))
                    else:
                        delta = 0.0
                    if delta < best_delta:
                        best_delta = delta
                        best_idx = idx
                t_s_idx = best_idx

        t_s = float(thresholds[t_s_idx])
        t_s2 = float(thresholds[t_s_idx + 1]) if t_s_idx + 1 < len(thresholds) else t_s
        return t_s, t_s2, thresholds.tolist()

    def manual_soma_setting(self, astrocyte_stack, initial_points=None, initial_types=None):
        self.ctx.log(self.ctx.t("log.soma_title"))

        z_dim, y_dim, x_dim = astrocyte_stack.shape
        astro = astrocyte_stack.astype(np.float32)
        astro_smooth = gaussian_filter(astro, sigma=0.5)
        nz = astro_smooth[astro_smooth > 0]
        if nz.size == 0:
            self.ctx.log(self.ctx.t("log.soma_no_points"))
            return np.zeros_like(astrocyte_stack, dtype=np.int32)
        p0, p99 = np.percentile(nz, [1, 99])
        astro_norm = np.clip((astro_smooth - p0) / (p99 - p0 + 1e-6), 0, 1)
        mip = np.max(astro_norm, axis=0)

        res = self.ctx.bridge._ask("soma_picker", mip=mip,
                                   astro_norm=astro_norm,
                                   initial_points=initial_points,
                                   initial_types=initial_types)

        if not res:
            self.ctx.log(self.ctx.t("log.soma_no_points"))
            return np.zeros_like(astrocyte_stack, dtype=np.int32)

        if isinstance(res, tuple) and len(res) == 2:
            points, point_types = res
        else:
            points, point_types = res, ["solid"] * len(res)

        self.ctx.log(self.ctx.t("log.soma_points_set", n=len(points)))

        click_centers = []
        click_types = []

        if self.ctx.scale_initialized:
            xy_um = float(self.ctx.xy_scale_um) * float(self.ctx.downsample)
        else:
            xy_um = 1.0

        z_window_px = max(3, int(5.0 / xy_um))

        for i, pt in enumerate(points):
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
                    max(0, soma_py - z_window_px):min(y_dim, soma_py + z_window_px + 1),
                    max(0, soma_px - z_window_px):min(x_dim, soma_px + z_window_px + 1)]
                    if window.size > 0:
                        score = np.median(window)
                        if score > best_score_click:
                            best_score_click = score
                            best_z_click = z

            click_centers.append((best_z_click, soma_py, soma_px))
            click_types.append(point_types[i] if i < len(point_types) else "solid")

        while True:
            res_confirm = self.ctx.bridge._ask(
                "soma_confirm",
                astro_norm=astro_norm,
                click_centers=click_centers,
                click_types=click_types,
            )

            if res_confirm == "to_manual":
                initial_pts = [(cy, cx, cz) for (cz, cy, cx) in click_centers]
                manual_res = self.ctx.bridge._ask(
                    "soma_picker",
                    mip=mip,
                    astro_norm=astro_norm,
                    initial_points=initial_pts,
                    initial_types=click_types,
                )
                if not manual_res:
                    self.ctx.log(self.ctx.t("log.soma_cancelled"))
                    return np.zeros_like(astrocyte_stack, dtype=np.int32)

                if isinstance(manual_res, tuple) and len(manual_res) == 2:
                    points, point_types = manual_res
                else:
                    points, point_types = manual_res, ["solid"] * len(manual_res)

                click_centers = []
                click_types = []
                for i, pt in enumerate(points):
                    py, px, pz = pt
                    py = max(0, min(y_dim - 1, py))
                    px = max(0, min(x_dim - 1, px))
                    pz = max(0, min(z_dim - 1, pz))
                    click_centers.append((pz, py, px))
                    click_types.append(point_types[i] if i < len(point_types) else "solid")
                continue

            if not res_confirm:
                self.ctx.log(self.ctx.t("log.soma_cancelled"))
                return np.zeros_like(astrocyte_stack, dtype=np.int32)

            break

        self.ctx.log(self.ctx.t("log.soma_starting"))

        import heapq
        from scipy.ndimage import label, binary_dilation, distance_transform_edt, binary_fill_holes
        from skimage.filters import threshold_otsu
        from skimage.segmentation import watershed
        from scipy.ndimage import label as nd_label
        from scipy.ndimage import distance_transform_edt as nd_edt
        from skimage.segmentation import watershed as sk_watershed

        grain_map = getattr(self.ctx, "_grain_index_map", None)
        soma_segmented_3d = np.zeros_like(astrocyte_stack, dtype=np.int32)

        if self.ctx.scale_initialized:
            xy_um = float(self.ctx.xy_scale_um) * float(self.ctx.downsample)
            z_um = float(self.ctx.z_scale_um)
            if z_um <= 0:
                z_um = 1.0
        else:
            xy_um = 1.0
            z_um = 1.0

        window_um = 30.0
        half_win_xy = int(round((window_um * 0.5) / max(xy_um, 1e-6)))
        half_win_xy = max(4, min(half_win_xy, max(y_dim, x_dim)))

        max_radius_xy = min(int(round(15.0 / max(xy_um, 1e-6))), 30)
        max_radius_xy = max(2, min(max_radius_xy, half_win_xy))

        max_center_dist_px = max(1, int(round(5.0 / max(xy_um, 1e-6))))
        max_center_dist_px = min(max_center_dist_px, half_win_xy)

        max_r_soma_um = 8.0
        max_r_soma_px = max(2, int(round(max_r_soma_um / max(xy_um, 1e-6))))
        max_r_soma_px = min(max_r_soma_px, max_radius_xy, half_win_xy)

        circle_mask_cache = {}
        for rr in range(1, max_r_soma_px + 1):
            yy_r, xx_r = np.ogrid[-rr:rr+1, -rr:rr+1]
            circle_mask_cache[rr] = (yy_r**2 + xx_r**2) <= rr**2

        all_candidate_masks = {}
        all_point_markers = []
        handled_donuts = []

        for i, (cz, cy, cx) in enumerate(click_centers, 1):
            if not (0 <= cz < z_dim and 0 <= cy < y_dim and 0 <= cx < x_dim):
                continue

            stype = click_types[i - 1] if i - 1 < len(click_types) else "solid"

            if stype == "donut":
                max_radius_z_d = 2
                max_radius_xy_d = 3
                donut_mask = np.zeros_like(astrocyte_stack, dtype=bool)
                for dz in range(-max_radius_z_d, max_radius_z_d + 1):
                    nz_ = cz + dz
                    if not (0 <= nz_ < z_dim):
                        continue
                    for dy in range(-max_radius_xy_d, max_radius_xy_d + 1):
                        ny = cy + dy
                        if not (0 <= ny < y_dim):
                            continue
                        for dx in range(-max_radius_xy_d, max_radius_xy_d + 1):
                            nx = cx + dx
                            if not (0 <= nx < x_dim):
                                continue
                            if (dz/max_radius_z_d)**2 + (dy/max_radius_xy_d)**2 + (dx/max_radius_xy_d)**2 <= 1.0:
                                donut_mask[nz_, ny, nx] = True
                handled_donuts.append((i, donut_mask))
                continue

            y0 = max(0, cy - half_win_xy)
            y1 = min(y_dim, cy + half_win_xy + 1)
            x0 = max(0, cx - half_win_xy)
            x1 = min(x_dim, cx + half_win_xy + 1)

            vol = astro_norm[:, y0:y1, x0:x1]
            h_local, w_local = y1 - y0, x1 - x0
            seed_y_local = cy - y0
            seed_x_local = cx - x0

            slice_masks = []

            r_search = max(1, int(round(max_center_dist_px)))
            yy_grid, xx_grid = np.mgrid[-r_search:r_search + 1, -r_search:r_search + 1]
            dist2_grid = yy_grid * yy_grid + xx_grid * xx_grid
            offsets_raw = np.argwhere(dist2_grid <= (r_search * r_search))
            offsets = sorted(offsets_raw.tolist(), key=lambda t: (t[0] - r_search) ** 2 + (t[1] - r_search) ** 2)

            for z_idx in range(z_dim):
                slice_data = vol[z_idx]
                if slice_data.size == 0 or np.ptp(slice_data) < 0.01:
                    slice_masks.append(np.zeros((h_local, w_local), dtype=bool))
                    continue

                nz_slice = slice_data[slice_data > 0]
                if nz_slice.size < 20:
                    slice_masks.append(np.zeros((h_local, w_local), dtype=bool))
                    continue

                try:
                    thr = float(threshold_otsu(nz_slice))
                except Exception:
                    thr = float(np.percentile(nz_slice, 50))

                binary_slice = (slice_data >= thr)

                ys, xs = np.where(binary_slice)
                if len(ys) == 0:
                    slice_masks.append(np.zeros((h_local, w_local), dtype=bool))
                    continue

                seed_y_c = max(0, min(h_local - 1, seed_y_local))
                seed_x_c = max(0, min(w_local - 1, seed_x_local))

                best_r_here = 0

                for r in range(max_r_soma_px, 0, -1):
                    found_for_r = False
                    cmask = circle_mask_cache[r]
                    for (oy, ox) in offsets:
                        py_c = seed_y_c + (oy - r_search)
                        px_c = seed_x_c + (ox - r_search)
                        if py_c < r or py_c >= h_local - r or px_c < r or px_c >= w_local - r:
                            continue
                        dy_c = py_c - seed_y_local
                        dx_c = px_c - seed_x_local
                        if (dy_c * dy_c + dx_c * dx_c) > (max_center_dist_px * max_center_dist_px):
                            continue

                        patch = binary_slice[py_c-r:py_c+r+1, px_c-r:px_c+r+1]
                        if np.all(patch[cmask]):
                            best_r_here = r
                            found_for_r = True
                            break
                    if found_for_r:
                        break

                marker_vol_slice = np.zeros((h_local, w_local), dtype=bool)

                if best_r_here > 0:
                    min_r = max(1, int(best_r_here * 0.75))
                    for r in range(best_r_here, min_r - 1, -1):
                        cmask = circle_mask_cache[r]
                        for (oy, ox) in offsets:
                            py_c = seed_y_c + (oy - r_search)
                            px_c = seed_x_c + (ox - r_search)
                            if py_c < r or py_c >= h_local - r or px_c < r or px_c >= w_local - r:
                                continue
                            dy_c = py_c - seed_y_local
                            dx_c = px_c - seed_x_local
                            if (dy_c * dy_c + dx_c * dx_c) > (max_center_dist_px * max_center_dist_px):
                                continue

                            patch = binary_slice[py_c-r:py_c+r+1, px_c-r:px_c+r+1]
                            if np.all(patch[cmask]):
                                marker_vol_slice[py_c-r:py_c+r+1, px_c-r:px_c+r+1][cmask] = True

                slice_masks.append(marker_vol_slice)

            marker_vol = np.stack(slice_masks, axis=0)
            full_mask_candidate = np.zeros_like(astrocyte_stack, dtype=bool)
            full_mask_candidate[:, y0:y1, x0:x1] = marker_vol

            if np.any(full_mask_candidate):
                all_candidate_masks[i] = full_mask_candidate
                all_point_markers.append((i, cz, cy, cx))

        for i, donut_mask in handled_donuts:
            soma_segmented_3d[donut_mask] = i

        if all_candidate_masks:
            union = np.zeros_like(astrocyte_stack, dtype=bool)
            for m in all_candidate_masks.values():
                union |= m

            markers_img = np.zeros_like(astrocyte_stack, dtype=np.int32)
            for i, cz, cy, cx in all_point_markers:
                if 0 <= cz < z_dim and 0 <= cy < y_dim and 0 <= cx < x_dim:
                    markers_img[cz, cy, cx] = i

            if np.any(markers_img > 0):
                sampling = (z_um, xy_um, xy_um)
                dist = nd_edt(~union, sampling=sampling).astype(np.float32)
                labeled_ws = sk_watershed(-dist, markers_img, mask=union)
                for i in all_candidate_masks.keys():
                    soma_segmented_3d[labeled_ws == i] = i
            else:
                for i, m in all_candidate_masks.items():
                    soma_segmented_3d[m & (soma_segmented_3d == 0)] = i

        self.ctx._soma_labels = soma_segmented_3d
        self.ctx._last_click_centers = click_centers
        self.ctx._last_click_types = click_types
        return soma_segmented_3d.astype(np.int32)

    def auto_soma_setting(self, astrocyte_stack, astro_norm=None, initial_points=None, initial_types=None):
        self.ctx.log(self.ctx.t("log.soma_auto_starting"))

        z_dim, y_dim, x_dim = astrocyte_stack.shape

        if astro_norm is None:
            astro = astrocyte_stack.astype(np.float32)
            astro_smooth = gaussian_filter(astro, sigma=0.5)
            nz = astro_smooth[astro_smooth > 0]
            if nz.size < 100:
                self.ctx.log(self.ctx.t("log.soma_auto_none"))
                return None
            p0, p99 = np.percentile(nz, [1, 99])
            astro_norm = np.clip((astro_smooth - p0) / (p99 - p0 + 1e-6), 0, 1)

        mip = np.max(astro_norm, axis=0)

        points = []
        types = []
        if initial_points:
            z_max_allowed = z_dim - 1
            for i, pt in enumerate(initial_points):
                if isinstance(pt, (tuple, list)) and len(pt) == 3:
                    py_, px_, pz_ = pt
                    py_ = max(0, min(y_dim - 1, int(py_)))
                    px_ = max(0, min(x_dim - 1, int(px_)))
                    pz_ = max(0, min(z_max_allowed, int(pz_)))
                    points.append((py_, px_, pz_))
                    types.append(initial_types[i] if initial_types and i < len(initial_types) else "solid")

        if not points:
            t_s_pair = self._multilevel_soma_threshold(mip)
            if t_s_pair is None:
                self.ctx.log(self.ctx.t("log.soma_auto_none"))
                return None
            t_s, t_s2, all_thresholds = t_s_pair

            weight_mip = None
            weight_map = getattr(self.ctx, "_weight_map", None)
            if weight_map is not None and weight_map.shape == astro_norm.shape:
                weight_mip = np.max(weight_map, axis=0)

            peak_xy, peak_scores = self._soma_peaks_2d(
                mip, astro_norm, all_thresholds, weight_mip=weight_mip
            )

            if not peak_xy:
                self.ctx.log(self.ctx.t("log.soma_auto_none"))
                return None

            candidates = []
            for (py, px), sc in zip(peak_xy, peak_scores):
                z_layers = self._z_layers_for_soma(astro_norm, py, px)
                if not z_layers:
                    z_layers = [int(np.argmax(astro_norm[:, py, px]))]
                for pz in z_layers:
                    candidates.append({
                        "pz": int(pz), "py": int(py), "px": int(px),
                        "score": float(sc),
                        "method": "existing_2d",
                    })

            if not candidates:
                self.ctx.log(self.ctx.t("log.soma_auto_none"))
                return None

            if self.ctx.scale_initialized:
                xy_um = float(self.ctx.xy_scale_um) * float(self.ctx.downsample)
                z_um = float(self.ctx.z_scale_um)
            else:
                xy_um = 1.0
                z_um = 1.0

            soma_diam_px_xy = 12.0 / max(xy_um, 1e-6)
            soma_diam_px_z = 12.0 / max(z_um, 1e-6)
            match_r_xy = max(2.5, soma_diam_px_xy * 0.7)
            match_r_z = max(1.5, soma_diam_px_z * 0.7)

            border_margin = max(5, int(10.0 / max(xy_um, 1e-6)))
            valid_peaks = []

            for c in candidates:
                if c["py"] < border_margin or c["py"] > (y_dim - border_margin) or \
                        c["px"] < border_margin or c["px"] > (x_dim - border_margin):
                    continue
                if c["pz"] <= 1 or c["pz"] >= (z_dim - 2):
                    continue
                valid_peaks.append(c)

            if not valid_peaks:
                self.ctx.log(self.ctx.t("log.soma_auto_none"))
                return None

            merged = []
            for p in valid_peaks:
                found = False
                for mp in merged:
                    dxy = math.hypot(mp["py"] - p["py"], mp["px"] - p["px"])
                    dz = abs(mp["pz"] - p["pz"])
                    if dxy <= match_r_xy and dz <= match_r_z:
                        if p["score"] > mp["score"]:
                            mp["score"] = p["score"]
                            mp["pz"] = p["pz"]
                            mp["py"] = p["py"]
                            mp["px"] = p["px"]
                        found = True
                        break
                if not found:
                    merged.append(p)

            merged.sort(key=lambda c: c["score"], reverse=True)

            centers = [(c["pz"], c["py"], c["px"]) for c in merged]
            scores = [c["score"] for c in merged]
            min_dist_um = 10.0
            kept_centers = self._nms_3d(centers, scores, min_dist_um)

            final_points = []
            for cen in kept_centers:
                pz, py, px = int(cen[0]), int(cen[1]), int(cen[2])
                final_points.append((py, px, pz))

            points = final_points
            types = ["solid"] * len(points)

        self.ctx.log(self.ctx.t("log.soma_auto_found", n=len(points)))

        while True:
            click_centers = [(pz, py, px) for (py, px, pz) in points]

            res = self.ctx.bridge._ask(
                "soma_confirm",
                astro_norm=astro_norm,
                click_centers=click_centers,
                click_types=types,
            )

            if res == "to_manual":
                manual_res = self.ctx.bridge._ask(
                    "soma_picker",
                    mip=mip,
                    astro_norm=astro_norm,
                    initial_points=points,
                    initial_types=types,
                )
                if not manual_res:
                    self.ctx.log(self.ctx.t("log.soma_cancelled"))
                    return None

                if isinstance(manual_res, tuple) and len(manual_res) == 2:
                    points, types = manual_res
                else:
                    points, types = manual_res, ["solid"] * len(manual_res)

                self.ctx.log(self.ctx.t("log.soma_points_set", n=len(points)))
                continue

            if not res:
                self.ctx.log(self.ctx.t("log.soma_cancelled"))
                return None

            break

        self.ctx.log(self.ctx.t("log.soma_starting"))

        click_centers = []
        click_types = []
        for i, pt in enumerate(points, 1):
            py, px, pz = pt
            py = max(0, min(y_dim - 1, py))
            px = max(0, min(x_dim - 1, px))
            pz = max(0, min(z_dim - 1, pz))
            click_centers.append((pz, py, px))
            click_types.append(types[i - 1] if types and i - 1 < len(types) else "solid")

        from scipy.ndimage import label, binary_dilation, distance_transform_edt
        from skimage.filters import threshold_otsu
        from skimage.segmentation import watershed
        from scipy.ndimage import label as nd_label
        from scipy.ndimage import distance_transform_edt as nd_edt
        from skimage.segmentation import watershed as sk_watershed

        grain_map = getattr(self.ctx, "_grain_index_map", None)
        soma_segmented_3d = np.zeros_like(astrocyte_stack, dtype=np.int32)

        if self.ctx.scale_initialized:
            xy_um = float(self.ctx.xy_scale_um) * float(self.ctx.downsample)
            z_um = float(self.ctx.z_scale_um)
            if z_um <= 0:
                z_um = 1.0
        else:
            xy_um = 1.0
            z_um = 1.0

        window_um = 30.0
        half_win_xy = int(round((window_um * 0.5) / max(xy_um, 1e-6)))
        half_win_xy = max(4, min(half_win_xy, max(y_dim, x_dim)))

        max_radius_xy = min(int(round(15.0 / max(xy_um, 1e-6))), 30)
        max_radius_xy = max(2, min(max_radius_xy, half_win_xy))

        max_center_dist_px = max(1, int(round(5.0 / max(xy_um, 1e-6))))
        max_center_dist_px = min(max_center_dist_px, half_win_xy)

        max_r_soma_um = 8.0
        max_r_soma_px = max(2, int(round(max_r_soma_um / max(xy_um, 1e-6))))
        max_r_soma_px = min(max_r_soma_px, max_radius_xy, half_win_xy)

        circle_mask_cache = {}
        for rr in range(1, max_r_soma_px + 1):
            yy_r, xx_r = np.ogrid[-rr:rr+1, -rr:rr+1]
            circle_mask_cache[rr] = (yy_r**2 + xx_r**2) <= rr**2

        all_candidate_masks = {}
        all_point_markers = []
        handled_donuts = []

        for i, (cz, cy, cx) in enumerate(click_centers, 1):
            if not (0 <= cz < z_dim and 0 <= cy < y_dim and 0 <= cx < x_dim):
                continue

            stype = click_types[i - 1] if i - 1 < len(click_types) else "solid"

            if stype == "donut":
                max_radius_z_d = 2
                max_radius_xy_d = 3
                donut_mask = np.zeros_like(astrocyte_stack, dtype=bool)
                for dz in range(-max_radius_z_d, max_radius_z_d + 1):
                    nz_ = cz + dz
                    if not (0 <= nz_ < z_dim):
                        continue
                    for dy in range(-max_radius_xy_d, max_radius_xy_d + 1):
                        ny = cy + dy
                        if not (0 <= ny < y_dim):
                            continue
                        for dx in range(-max_radius_xy_d, max_radius_xy_d + 1):
                            nx = cx + dx
                            if not (0 <= nx < x_dim):
                                continue
                            if (dz/max_radius_z_d)**2 + (dy/max_radius_xy_d)**2 + (dx/max_radius_xy_d)**2 <= 1.0:
                                donut_mask[nz_, ny, nx] = True
                handled_donuts.append((i, donut_mask))
                continue

            y0 = max(0, cy - half_win_xy)
            y1 = min(y_dim, cy + half_win_xy + 1)
            x0 = max(0, cx - half_win_xy)
            x1 = min(x_dim, cx + half_win_xy + 1)

            vol = astro_norm[:, y0:y1, x0:x1]
            h_local, w_local = y1 - y0, x1 - x0
            seed_y_local = cy - y0
            seed_x_local = cx - x0

            slice_masks = []

            r_search = max(1, int(round(max_center_dist_px)))
            yy_grid, xx_grid = np.mgrid[-r_search:r_search + 1, -r_search:r_search + 1]
            dist2_grid = yy_grid * yy_grid + xx_grid * xx_grid
            offsets_raw = np.argwhere(dist2_grid <= (r_search * r_search))
            offsets = sorted(offsets_raw.tolist(), key=lambda t: (t[0] - r_search) ** 2 + (t[1] - r_search) ** 2)

            for z_idx in range(z_dim):
                slice_data = vol[z_idx]
                if slice_data.size == 0 or np.ptp(slice_data) < 0.01:
                    slice_masks.append(np.zeros((h_local, w_local), dtype=bool))
                    continue

                nz_slice = slice_data[slice_data > 0]
                if nz_slice.size < 20:
                    slice_masks.append(np.zeros((h_local, w_local), dtype=bool))
                    continue

                try:
                    thr = float(threshold_otsu(nz_slice))
                except Exception:
                    thr = float(np.percentile(nz_slice, 50))

                binary_slice = (slice_data >= thr)

                ys, xs = np.where(binary_slice)
                if len(ys) == 0:
                    slice_masks.append(np.zeros((h_local, w_local), dtype=bool))
                    continue

                seed_y_c = max(0, min(h_local - 1, seed_y_local))
                seed_x_c = max(0, min(w_local - 1, seed_x_local))

                best_r_here = 0

                for r in range(max_r_soma_px, 0, -1):
                    found_for_r = False
                    cmask = circle_mask_cache[r]
                    for (oy, ox) in offsets:
                        py_c = seed_y_c + (oy - r_search)
                        px_c = seed_x_c + (ox - r_search)
                        if py_c < r or py_c >= h_local - r or px_c < r or px_c >= w_local - r:
                            continue
                        dy_c = py_c - seed_y_local
                        dx_c = px_c - seed_x_local
                        if (dy_c * dy_c + dx_c * dx_c) > (max_center_dist_px * max_center_dist_px):
                            continue

                        patch = binary_slice[py_c-r:py_c+r+1, px_c-r:px_c+r+1]
                        if np.all(patch[cmask]):
                            best_r_here = r
                            found_for_r = True
                            break
                    if found_for_r:
                        break

                marker_vol_slice = np.zeros((h_local, w_local), dtype=bool)

                if best_r_here > 0:
                    min_r = max(1, int(best_r_here * 0.75))
                    for r in range(best_r_here, min_r - 1, -1):
                        cmask = circle_mask_cache[r]
                        for (oy, ox) in offsets:
                            py_c = seed_y_c + (oy - r_search)
                            px_c = seed_x_c + (ox - r_search)
                            if py_c < r or py_c >= h_local - r or px_c < r or px_c >= w_local - r:
                                continue
                            dy_c = py_c - seed_y_local
                            dx_c = px_c - seed_x_local
                            if (dy_c * dy_c + dx_c * dx_c) > (max_center_dist_px * max_center_dist_px):
                                continue

                            patch = binary_slice[py_c-r:py_c+r+1, px_c-r:px_c+r+1]
                            if np.all(patch[cmask]):
                                marker_vol_slice[py_c-r:py_c+r+1, px_c-r:px_c+r+1][cmask] = True

                slice_masks.append(marker_vol_slice)

            marker_vol = np.stack(slice_masks, axis=0)
            full_mask_candidate = np.zeros_like(astrocyte_stack, dtype=bool)
            full_mask_candidate[:, y0:y1, x0:x1] = marker_vol

            if np.any(full_mask_candidate):
                all_candidate_masks[i] = full_mask_candidate
                all_point_markers.append((i, cz, cy, cx))

        for i, donut_mask in handled_donuts:
            soma_segmented_3d[donut_mask] = i

        if all_candidate_masks:
            union = np.zeros_like(astrocyte_stack, dtype=bool)
            for m in all_candidate_masks.values():
                union |= m

            markers_img = np.zeros_like(astrocyte_stack, dtype=np.int32)
            for i, cz, cy, cx in all_point_markers:
                if 0 <= cz < z_dim and 0 <= cy < y_dim and 0 <= cx < x_dim:
                    markers_img[cz, cy, cx] = i

            if np.any(markers_img > 0):
                sampling = (z_um, xy_um, xy_um)
                dist = nd_edt(~union, sampling=sampling).astype(np.float32)
                labeled_ws = sk_watershed(-dist, markers_img, mask=union)
                for i in all_candidate_masks.keys():
                    soma_segmented_3d[labeled_ws == i] = i
            else:
                for i, m in all_candidate_masks.items():
                    soma_segmented_3d[m & (soma_segmented_3d == 0)] = i

        self.ctx._soma_labels = soma_segmented_3d
        self.ctx._last_click_centers = click_centers
        self.ctx._last_click_types = click_types
        return soma_segmented_3d.astype(np.int32)

    def _soma_peaks_2d(self, mip, astro_norm, thresholds, weight_mip=None):
        if not thresholds:
            return [], []

        h, w = mip.shape
        z_dim = astro_norm.shape[0]
        is_3d = (z_dim > 1)

        mip_norm = mip.astype(np.float32)
        nz = mip_norm[mip_norm > 0]
        if nz.size > 0:
            p2, p98 = np.percentile(nz, [2, 98])
            mip_norm = np.clip((mip_norm - p2) / (p98 - p2 + 1e-6), 0, 1)

        sigma_bg = max(15, min(h, w) // 15)
        mip_bg = uniform_filter(mip_norm, size=sigma_bg, mode='reflect')
        mip_bg = np.maximum(mip_bg, 1e-6)
        local_contrast = mip_norm / mip_bg

        weight_norm = None
        if weight_mip is not None and weight_mip.shape == (h, w):
            wmax = float(np.max(weight_mip))
            if wmax > 0:
                weight_norm = (weight_mip / wmax).astype(np.float32)

        def _find_peaks_at(t):
            mask = (mip >= float(t))
            if not np.any(mask):
                return []
            mask = remove_small_objects(mask, min_size=30)
            if not np.any(mask):
                return []
            edt = distance_transform_edt(mask).astype(np.float32)
            if edt.max() <= 0:
                return []
            coords = np.argwhere(
                (edt > 0) &
                (edt == maximum_filter(edt, size=9))
            )
            peaks = []
            for (py, px) in coords:
                r_core = float(edt[py, px])
                if r_core < 1.0:
                    continue
                peaks.append((int(py), int(px), r_core))
            return peaks

        per_level_peaks = []
        for t in thresholds:
            per_level_peaks.append(_find_peaks_at(t))

        base_level_idx = len(thresholds) - 1
        for i in range(len(thresholds) - 1, -1, -1):
            if len(per_level_peaks[i]) > 0:
                base_level_idx = i
                break

        base_peaks = per_level_peaks[base_level_idx]
        if not base_peaks:
            base_level_idx = 0
            base_peaks = per_level_peaks[0]
        if not base_peaks:
            return [], []

        base_t = float(thresholds[base_level_idx])
        base_mask = (mip >= base_t)
        base_mask = remove_small_objects(base_mask, min_size=30)

        def _isotropy_at_scale(py, px, scale_radius):
            r = max(2, int(round(scale_radius)))
            y0, y1 = max(0, py - r), min(h, py + r + 1)
            x0, x1 = max(0, px - r), min(w, px + r + 1)
            win = base_mask[y0:y1, x0:x1]
            if not np.any(win):
                return 0.0
            center_y = py - y0
            center_x = px - x0
            directions = [(0,1),(1,1),(1,0),(1,-1),(0,-1),(-1,-1),(-1,0),(-1,1)]
            decay_lengths = []
            for (dy, dx) in directions:
                length = 0.0
                step = 1.0
                while step <= r:
                    yy = int(round(center_y + dy * step))
                    xx = int(round(center_x + dx * step))
                    if not (0 <= yy < win.shape[0] and 0 <= xx < win.shape[1]):
                        break
                    if not win[yy, xx]:
                        break
                    length = step
                    step += 1.0
                decay_lengths.append(length)
            valid = [L for L in decay_lengths if L > 0]
            if len(valid) < 4:
                return 0.0
            min_decay = min(valid)
            max_decay = max(valid)
            if max_decay <= 0:
                return 0.0
            return float(min_decay / max_decay)

        def _volume_score(py, px, r_core):
            r = max(2, int(round(r_core * 1.5)))
            y0, y1 = max(0, py - r), min(h, py + r + 1)
            x0, x1 = max(0, px - r), min(w, px + r + 1)

            z_profile = astro_norm[:, y0:y1, x0:x1].mean(axis=(1, 2))
            if z_profile.size == 0:
                return 0.0
            prof_n = z_profile / (z_profile.max() + 1e-6)
            above = prof_n >= 0.5
            if not np.any(above):
                return 0.0

            z_indices = np.where(above)[0]
            z_span = int(z_indices.max() - z_indices.min() + 1)
            xy_span = min(y1 - y0, x1 - x0)

            if z_span < 2:
                return 0.0
            if z_span > xy_span * 3:
                return 0.0

            core_mask = astro_norm[z_indices.min():z_indices.max()+1, y0:y1, x0:x1]
            if core_mask.size == 0:
                return 0.0
            high_frac = float(np.mean(core_mask >= 0.5 * core_mask.max()))

            return high_frac

        raw_candidates = []
        for (py, px, r_core) in base_peaks:
            iso1 = _isotropy_at_scale(py, px, r_core * 1.0)
            iso2 = _isotropy_at_scale(py, px, r_core * 2.0)
            iso3 = _isotropy_at_scale(py, px, r_core * 3.0)
            min_iso = min(iso1, iso2, iso3) if (iso1 > 0 and iso2 > 0 and iso3 > 0) else 0.0

            persistence = 1
            for step in range(1, len(thresholds)):
                idx = base_level_idx - step
                if idx < 0:
                    break
                other = per_level_peaks[idx]
                if not other:
                    break
                tol = max(2.0, r_core)
                found = False
                for (oy, ox, orr) in other:
                    if math.hypot(py - oy, px - ox) <= tol:
                        found = True
                        break
                if found:
                    persistence += 1
                else:
                    break

            contrast_val = float(local_contrast[py, px])
            weight_val = float(weight_norm[py, px]) if weight_norm is not None else 1.0
            vol_score = _volume_score(py, px, r_core)

            raw_candidates.append({
                'py': int(py),
                'px': int(px),
                'r_core': r_core,
                'min_iso': min_iso,
                'persistence': persistence,
                'contrast': contrast_val,
                'weight': weight_val,
                'vol_score': vol_score,
            })

        if not raw_candidates:
            return [], []

        r_vals = np.array([c['r_core'] for c in raw_candidates], dtype=np.float32)
        iso_vals = np.array([c['min_iso'] for c in raw_candidates], dtype=np.float32)
        pers_vals = np.array([c['persistence'] for c in raw_candidates], dtype=np.float32)
        contrast_vals = np.array([c['contrast'] for c in raw_candidates], dtype=np.float32)
        weight_vals = np.array([c['weight'] for c in raw_candidates], dtype=np.float32)
        vol_vals = np.array([c['vol_score'] for c in raw_candidates], dtype=np.float32)

        r_ref = float(np.median(r_vals))
        iso_ref = float(np.median(iso_vals))
        pers_ref = float(np.median(pers_vals))
        contrast_ref = float(np.median(contrast_vals))
        weight_ref = float(np.median(weight_vals))
        vol_ref = float(np.median(vol_vals[vol_vals > 0])) if np.any(vol_vals > 0) else 0.0

        def _soft_score(c):
            s = 0
            if c['r_core'] >= r_ref * 0.5:
                s += 1
            if c['min_iso'] >= iso_ref * 0.5:
                s += 1
            if c['persistence'] >= max(1, pers_ref - 1):
                s += 1
            if c['contrast'] >= contrast_ref * 0.5:
                s += 1
            if c['weight'] >= weight_ref * 0.5:
                s += 1
            if c['vol_score'] >= vol_ref * 0.5:
                s += 1
            return s

        kept_soft = [c for c in raw_candidates if _soft_score(c) >= 2]

        if not kept_soft:
            kept_soft = sorted(raw_candidates,
                               key=lambda c: c['vol_score'] * c['r_core'],
                               reverse=True)[:5]

        if vol_ref > 0:
            kept_hard = [c for c in kept_soft if c['vol_score'] >= vol_ref * 0.3]
            if kept_hard:
                kept_soft = kept_hard

        kept_soft.sort(key=lambda c: c['vol_score'] * c['r_core'] * c['weight'], reverse=True)

        ys = [c['py'] for c in kept_soft]
        xs = [c['px'] for c in kept_soft]
        scores = [c['vol_score'] * c['r_core'] * c['weight'] for c in kept_soft]
        return list(zip(ys, xs)), scores

    def _z_layers_for_soma(self, astro_norm, py, px):
        z_dim, y_dim, x_dim = astro_norm.shape
        r = 3
        y0, y1 = max(0, py - r), min(y_dim, py + r + 1)
        x0, x1 = max(0, px - r), min(x_dim, px + r + 1)
        profile = astro_norm[:, y0:y1, x0:x1].mean(axis=(1, 2))
        if profile.size == 0 or profile.max() <= 0:
            return []
        prof_n = profile / (profile.max() + 1e-6)

        prominence = 0.15
        peaks, _ = find_peaks(prof_n, prominence=prominence)
        if len(peaks) == 0:
            peaks = np.array([int(np.argmax(prof_n))])

        if self.ctx.scale_initialized:
            z_um = float(self.ctx.z_scale_um)
        else:
            z_um = 1.0

        min_sep_layers = max(1, int(round(2.5 / max(z_um, 1e-6))))
        kept = []
        for p in peaks:
            p = int(p)
            if not kept:
                kept.append(p)
                continue
            if p - kept[-1] >= min_sep_layers:
                kept.append(p)
            else:
                if prof_n[p] > prof_n[kept[-1]]:
                    kept[-1] = p
        return [int(z) for z in kept]

    def _nms_3d(self, centers, scores, min_dist_px_z):
        if not centers:
            return []
        order = np.argsort(scores)[::-1]
        z_um = float(self.ctx.z_scale_um) if self.ctx.scale_initialized else 1.0
        xy_um = float(self.ctx.xy_scale_um) * float(self.ctx.downsample) if self.ctx.scale_initialized else 1.0
        keep = []
        for i in order:
            ci = centers[i]
            ok = True
            for j in keep:
                cj = centers[j]
                dz = (ci[0] - cj[0]) * z_um
                dy = (ci[1] - cj[1]) * xy_um
                dx = (ci[2] - cj[2]) * xy_um
                if math.sqrt(dz*dz + dy*dy + dx*dx) < min_dist_px_z:
                    ok = False
                    break
            if ok:
                keep.append(i)
        return [centers[i] for i in keep]

    def _detect_donut_and_fill(self, astro_norm_3d, click_centers, click_types, soma_mask):
        from scipy.ndimage import distance_transform_edt
        from skimage.segmentation import watershed

        volume = np.asarray(astro_norm_3d, dtype=np.float32)
        z_dim, y_dim, x_dim = volume.shape

        if getattr(self.ctx, "has_dapi", False) and getattr(self.ctx, "nuclei_segmented", None) is not None:
            nuclei = self.ctx.nuclei_segmented
            unique_nucs = np.unique(nuclei[nuclei > 0])
            if len(unique_nucs) > 0:
                markers = np.zeros_like(nuclei, dtype=np.int32)
                for i, nid in enumerate(unique_nucs, 1):
                    markers[nuclei == nid] = i

                dist = distance_transform_edt(soma_mask)
                new_soma = watershed(-dist, markers, mask=(soma_mask > 0))

                final_soma = np.zeros_like(soma_mask, dtype=np.int32)
                for i, nid in enumerate(unique_nucs, 1):
                    orig_cid = int(np.unique(soma_mask[markers == i])[0])
                    final_soma[new_soma == i] = orig_cid

                soma_mask = final_soma

        xy_um = 1.0
        z_um = 1.0
        ctx = getattr(self, "ctx", None)
        if ctx is not None and getattr(ctx, "scale_initialized", False):
            try:
                xy_um = float(ctx.xy_scale_um) * float(ctx.downsample)
                z_um = float(ctx.z_scale_um)
            except (TypeError, ValueError, AttributeError):
                xy_um = 1.0
                z_um = 1.0
        if not np.isfinite(xy_um) or xy_um <= 0:
            xy_um = 1.0
        if not np.isfinite(z_um) or z_um <= 0:
            z_um = 1.0

        donut_indices = [i for i, t in enumerate(click_types) if t == "donut"]
        if not donut_indices:
            return soma_mask.copy()

        updated_soma = soma_mask.copy()
        raw_grown_masks = {}

        for idx in donut_indices:
            if idx >= len(click_centers):
                continue
            cz, cy, cx = click_centers[idx]
            cz = int(max(0, min(z_dim - 1, cz)))
            cy = int(max(0, min(y_dim - 1, cy)))
            cx = int(max(0, min(x_dim - 1, cx)))

            cid = int(updated_soma[cz, cy, cx])
            if cid <= 0:
                cid = idx + 1

            id_map = self._donut_tile_fill(
                volume, (cz, cy, cx), cid, xy_um,
                astro_norm_clean=astro_norm_3d
            )

            if id_map is None or not np.any(id_map > 0):
                continue

            mask_region = (id_map == cid)
            if not np.any(mask_region):
                continue

            raw_grown_masks[cid] = mask_region

        if not raw_grown_masks:
            return updated_soma

        sampling = (z_um, xy_um, xy_um)
        dist_maps = {}
        for cid, mask in raw_grown_masks.items():
            dist_maps[cid] = distance_transform_edt(~mask, sampling=sampling).astype(np.float32)

        all_grown_union = np.zeros_like(volume, dtype=bool)
        for mask in raw_grown_masks.values():
            all_grown_union |= mask

        assigned_cids = np.zeros_like(volume, dtype=np.int32)
        min_dists = np.full(volume.shape, np.inf, dtype=np.float32)

        for cid, dist_map in dist_maps.items():
            update_mask = all_grown_union & (dist_map < min_dists)
            assigned_cids[update_mask] = cid
            min_dists[update_mask] = dist_map[update_mask]

        for cid in raw_grown_masks.keys():
            updated_soma[assigned_cids == cid] = cid
            filled_cid_mask = binary_fill_holes(updated_soma == cid)
            updated_soma[filled_cid_mask] = cid

        return updated_soma

    def _donut_tile_fill(self, volume, center, cid, xy_um, astro_norm_clean=None):
        from skimage.filters import threshold_otsu
        from skimage.morphology import remove_small_objects
        from scipy.ndimage import (
            binary_dilation, binary_propagation, label, distance_transform_edt,
            gaussian_filter, maximum_filter, uniform_filter,
            binary_closing
        )

        volume = np.asarray(volume, dtype=np.float32)
        z_dim, y_dim, x_dim = volume.shape
        cz0, cy0, cx0 = map(int, center)

        if astro_norm_clean is None:
            astro_norm_clean = volume
        else:
            astro_norm_clean = np.asarray(astro_norm_clean, dtype=np.float32)

        xy_um = float(xy_um)
        if not np.isfinite(xy_um) or xy_um <= 0:
            xy_um = 1.0

        tile_um = 40.0
        margin_um = 5.0
        shift_um = max(10.0, 10.0 * xy_um)

        tile_px = max(16, int(round(tile_um / xy_um)))
        margin_px = max(3, int(round(margin_um / xy_um)))
        shift_px = max(4, int(round(shift_um / xy_um)))

        outer = tile_px + 2 * margin_px
        inner = tile_px

        if outer >= y_dim or outer >= x_dim:
            outer = max(8, min(y_dim, x_dim) - 2)
            margin_px = min(margin_px, max(1, (outer - 4) // 2))
            inner = max(2, outer - 2 * margin_px)

        import hashlib

        volume_digest = hashlib.blake2b(
            np.ascontiguousarray(volume).view(np.uint8),
            digest_size=16
        ).hexdigest()
        volume_hash = (tuple(volume.shape), str(volume.dtype), volume_digest)
        cache_key = (volume_hash, round(xy_um, 4), outer, inner)

        if getattr(self.ctx, "_donut_maps_cache", None) is None:
            self.ctx._donut_maps_cache = {}

        if cache_key in self.ctx._donut_maps_cache:
            hits_count_3d, donut_mask_bin, labeled_map, peak_coords = self.ctx._donut_maps_cache[cache_key]
        else:
            volume_mip = np.zeros_like(volume)
            for zz in range(z_dim):
                z_lo = max(0, zz - 1)
                z_hi = min(z_dim, zz + 2)
                volume_mip[zz] = np.max(volume[z_lo:z_hi], axis=0)

            ys = list(range(0, max(1, y_dim - outer + 1), shift_px))
            xs = list(range(0, max(1, x_dim - outer + 1), shift_px))
            if not ys or ys[-1] != y_dim - outer:
                ys.append(max(0, y_dim - outer))
            if not xs or xs[-1] != x_dim - outer:
                xs.append(max(0, x_dim - outer))

            hits_count_3d = np.zeros((z_dim, y_dim, x_dim), dtype=np.int32)

            for z in range(z_dim):
                sub = volume_mip[z]
                nz_all = sub[sub > 0]
                if nz_all.size < 50:
                    continue

                bg = uniform_filter(sub, size=max(15, outer // 2))
                sub_enhanced = sub / np.maximum(bg, 1e-6)

                for ty in ys:
                    for tx in xs:
                        tile = sub_enhanced[ty:ty + outer, tx:tx + outer]
                        tile_orig = sub[ty:ty + outer, tx:tx + outer]

                        valid_mask = tile_orig > 0
                        valid_holes = np.zeros(tile.shape, dtype=bool)

                        if np.count_nonzero(valid_mask) >= 30:
                            tile_smooth = gaussian_filter(tile, sigma=0.8)

                            nz_vals = tile_smooth[valid_mask]
                            if nz_vals.size >= 30 and np.ptp(nz_vals) > 0:
                                try:
                                    thr_otsu = float(threshold_otsu(nz_vals))
                                except Exception:
                                    thr_otsu = 1.0

                                thr = max(thr_otsu * 0.5, 0.1)
                                binary = (tile_smooth >= thr) & valid_mask

                                if np.any(binary):
                                    binary_d = binary_dilation(binary, iterations=1)
                                    background = ~binary_d
                                    border = np.zeros_like(background)
                                    border[0, :] = True
                                    border[-1, :] = True
                                    border[:, 0] = True
                                    border[:, -1] = True

                                    exterior = binary_propagation(border & background, mask=background)
                                    holes = background & ~exterior

                                    if np.any(holes):
                                        labeled_holes, n_holes = label(holes, structure=np.ones((3, 3)))

                                        for hid in range(1, n_holes + 1):
                                            candidate = labeled_holes == hid
                                            area = int(np.count_nonzero(candidate))

                                            if area < 15 or area > 0.4 * tile.size:
                                                continue

                                            py, px = np.where(candidate)
                                            if py.size == 0:
                                                continue

                                            h = int(py.max() - py.min() + 1)
                                            w = int(px.max() - px.min() + 1)
                                            if h < 3 or w < 3:
                                                continue

                                            box_area = h * w
                                            fill_ratio = area / max(box_area, 1)

                                            if fill_ratio < 0.25:
                                                continue

                                            valid_holes |= candidate

                                        valid_holes &= ~binary

                        if not np.any(valid_holes):
                            continue

                        y0 = margin_px
                        x0 = margin_px
                        y1 = min(y0 + inner, outer)
                        x1 = min(x0 + inner, outer)

                        holes_inner = valid_holes[y0:y1, x0:x1]
                        if not np.any(holes_inner):
                            continue

                        gy0, gx0 = ty + y0, tx + x0
                        gy1, gx1 = gy0 + holes_inner.shape[0], gx0 + holes_inner.shape[1]
                        hits_count_3d[z, gy0:gy1, gx0:gx1] += holes_inner.astype(np.int32)

            max_hits = int(np.max(hits_count_3d)) if hits_count_3d.size > 0 else 0
            if max_hits < 3:
                self.ctx._donut_maps_cache[cache_key] = (
                    hits_count_3d, np.zeros_like(hits_count_3d, dtype=bool), None, []
                )
                return None
            else:
                donut_mask_raw = hits_count_3d >= 3
                donut_mask_raw = remove_small_objects(donut_mask_raw, min_size=3)
                donut_mask_bin = binary_closing(
                    donut_mask_raw,
                    structure=np.ones((3, 3, 3), dtype=bool),
                    iterations=1,
                    border_value=0
                )
                donut_mask_bin |= donut_mask_raw
                donut_mask_bin = remove_small_objects(donut_mask_bin, min_size=3)

                labeled_map, _ = label(
                    donut_mask_bin,
                    structure=np.ones((3, 3, 3), dtype=bool)
                )

                peak_coords = []
                if np.any(donut_mask_bin):
                    slice_dist_2d_peak = max(3, int(round(10.0 / xy_um)))
                    for zz in range(z_dim):
                        sl = donut_mask_bin[zz]
                        if not np.any(sl):
                            continue
                        sl_dist = distance_transform_edt(sl)
                        sl_peaks = (sl_dist > 1.0) & (sl_dist == maximum_filter(sl_dist, size=slice_dist_2d_peak))
                        sl_peak_coords = np.argwhere(sl_peaks)
                        for (py, px) in sl_peak_coords:
                            peak_coords.append((zz, int(py), int(px)))
                    peak_coords = np.array(peak_coords) if peak_coords else np.zeros((0, 3), dtype=int)

                self.ctx._donut_maps_cache[cache_key] = (
                    hits_count_3d, donut_mask_bin, labeled_map, peak_coords
                )

        z_um = 1.0
        try:
            if getattr(self.ctx, "scale_initialized", False):
                z_um = float(self.ctx.z_scale_um)
        except (TypeError, ValueError, AttributeError):
            z_um = 1.0

        if not np.isfinite(z_um) or z_um <= 0:
            z_um = 1.0

        max_hits = int(np.max(hits_count_3d)) if hits_count_3d.size > 0 else 0
        if max_hits <= 0 or labeled_map is None or not np.any(donut_mask_bin):
            return None

        unique_ids = np.unique(labeled_map)
        unique_ids = unique_ids[unique_ids > 0]

        best_id = 0
        best_score = -float("inf")

        inside_cid = int(labeled_map[cz0, cy0, cx0])

        if inside_cid > 0:
            best_id = inside_cid
        else:
            for comp_id in unique_ids:
                coords = np.argwhere(labeled_map == comp_id)
                if not coords.size:
                    continue

                dz = (coords[:, 0] - cz0) * z_um
                dy = (coords[:, 1] - cy0) * xy_um
                dx = (coords[:, 2] - cx0) * xy_um
                dists_3d = np.sqrt(dz * dz + dy * dy + dx * dx)
                min_idx = int(np.argmin(dists_3d))
                dmin_3d = float(dists_3d[min_idx])

                z_slices = np.unique(coords[:, 0])
                n_z = int(len(z_slices))
                solidity_vals = []

                for zz in z_slices:
                    sl_mask = (labeled_map[zz] == comp_id)
                    props_2d = measure.regionprops(sl_mask.astype(np.uint8))
                    if props_2d:
                        solidity_vals.append(float(props_2d[0].solidity))

                mean_solidity = (
                    float(np.mean(solidity_vals)) if solidity_vals else 0.0
                )
                volume_um3 = float(coords.shape[0]) * (xy_um ** 2) * z_um

                score = (
                        -2.0 * dmin_3d
                        + 1.5 * n_z
                        + 2.0 * mean_solidity
                )
                accepted = dmin_3d <= 10.0

                if not accepted:
                    continue

                if score > best_score:
                    best_score = score
                    best_id = int(comp_id)

        if best_id == 0:
            return None

        soma_labels = getattr(self.ctx, "_soma_labels", None)
        soma_other_mask = np.zeros_like(volume, dtype=bool)
        if soma_labels is not None and soma_labels.shape == volume.shape:
            soma_other_mask = (soma_labels > 0) & (soma_labels != cid)

        seed = (labeled_map == best_id)

        grown_mask = self._grow_donut_masks_joint(
            volume, seed, cid, xy_um, z_um,
            soma_other_mask, astro_norm_clean=astro_norm_clean,
            center=(cz0, cy0, cx0)
        )

        if grown_mask is None or not np.any(grown_mask):
            return None

        id_map = np.zeros_like(volume, dtype=np.int32)
        id_map[grown_mask] = cid

        return id_map

    def _grow_donut_masks_joint(self, volume, seed_mask, cid, xy_um, z_um,
                                soma_other_mask, astro_norm_clean=None,
                                center=None):
        from scipy.ndimage import (
            distance_transform_edt, gaussian_filter, binary_fill_holes,
            binary_closing
        )
        from skimage.morphology import remove_small_objects
        from skimage.filters import threshold_otsu

        volume = np.asarray(volume, dtype=np.float32)
        z_dim, y_dim, x_dim = volume.shape

        if astro_norm_clean is None:
            astro_norm_clean = volume
        else:
            astro_norm_clean = np.asarray(astro_norm_clean, dtype=np.float32)

        seed_mask = np.asarray(seed_mask, dtype=bool)
        if not np.any(seed_mask):
            return None

        sampling = (z_um, xy_um, xy_um)
        step_um = 1.0
        max_iter = 60
        sigma_xy_local = max(3.0, 8.0 / max(xy_um, 1e-6))
        sigma_z_local = max(1.0, 4.0 / max(z_um, 1e-6))

        nz_vals = astro_norm_clean[astro_norm_clean > 0]
        if nz_vals.size < 50:
            return None

        try:
            background_thr = float(threshold_otsu(nz_vals))
        except Exception:
            background_thr = float(np.percentile(nz_vals, 50))

        background_mask = astro_norm_clean <= background_thr
        bg_mean = float(np.mean(astro_norm_clean[background_mask])) if np.any(background_mask) else 0.0
        bg_std = float(np.std(astro_norm_clean[background_mask])) if np.any(background_mask) else 0.0

        signal_thr = background_thr + 2.0 * bg_std

        seed_vals = astro_norm_clean[seed_mask]
        if seed_vals.size < 10:
            return None
        seed_p20 = float(np.percentile(seed_vals, 20))
        seed_p50 = float(np.percentile(seed_vals, 50))

        seed_thr = max(seed_p20, background_thr + 1.5 * bg_std)
        bright_thr = max(seed_p50, background_thr + 2.5 * bg_std)

        core_edt = distance_transform_edt(astro_norm_clean > bright_thr, sampling=sampling)
        seed_edt = distance_transform_edt(seed_mask, sampling=sampling)
        r_seed = float(np.percentile(seed_edt[seed_mask], 90)) if np.any(seed_mask) else 0.0
        if not np.isfinite(r_seed) or r_seed <= 0:
            r_seed = max(2.0, 5.0 / max(xy_um, 1e-6))

        thin_max_um = max(0.75, 0.4 * r_seed)
        thin = core_edt < thin_max_um

        passable_iso = (astro_norm_clean >= seed_thr)
        passable_iso &= (~soma_other_mask)
        passable_iso &= (~thin)

        if center is not None:
            cz0, cy0, cx0 = center
            max_radius_um = 3.0 * r_seed
            max_radius_z_um = 2.0 * r_seed

            zz, yy, xx = np.ogrid[:z_dim, :y_dim, :x_dim]
            dz = (zz - cz0) * z_um
            dy = (yy - cy0) * xy_um
            dx = (xx - cx0) * xy_um
            ellipsoid = (dz ** 2 / max(max_radius_z_um, 1e-6) ** 2 +
                         dy ** 2 / max(max_radius_um, 1e-6) ** 2 +
                         dx ** 2 / max(max_radius_um, 1e-6) ** 2) <= 1.0
            passable_iso &= ellipsoid

        anchor_dist = distance_transform_edt(~seed_mask, sampling=sampling).astype(np.float32)
        anchor_um = max(2.0, 1.5 * r_seed)

        volume_smooth = gaussian_filter(
            astro_norm_clean,
            sigma=(sigma_z_local, sigma_xy_local, sigma_xy_local),
            mode='reflect'
        )

        current = seed_mask.copy()
        med_current = float(np.median(astro_norm_clean[current]))
        local_floor = max(seed_p20, med_current * 0.7)

        for _ in range(max_iter):
            if not np.any(current):
                break

            d_current = distance_transform_edt(~current, sampling=sampling).astype(np.float32)
            d_current = gaussian_filter(d_current, sigma=(0.5, 0.5, 0.5))

            grow_iso = passable_iso & (d_current <= step_um) & (d_current > 0)
            iso_cand = grow_iso & (~current)

            aniso_cand = np.zeros_like(current)
            near_anchor = (anchor_dist <= anchor_um) & (anchor_dist > 0)
            if np.any(near_anchor):
                local_med = float(np.median(volume_smooth[near_anchor & (astro_norm_clean > 0)]))
                aniso_thr = max(local_floor, local_med * 0.75)
                bright_near = near_anchor & (astro_norm_clean >= aniso_thr)
                aniso_cand = bright_near & (~current) & (~soma_other_mask)

                if np.any(aniso_cand):
                    temp_aniso = current | aniso_cand
                    cand_vals = astro_norm_clean[aniso_cand]
                    cur_vals = astro_norm_clean[current]

                    if np.any(cand_vals) and np.any(cur_vals):
                        p50_cand = float(np.percentile(cand_vals, 50))
                        p15_cur = float(np.percentile(cur_vals, 15))
                        p10_temp = float(np.percentile(astro_norm_clean[temp_aniso], 10))
                        p10_cur = float(np.percentile(cur_vals, 10))

                        if (p50_cand < p15_cur) or (p10_temp < p10_cur * 0.95):
                            aniso_cand = np.zeros_like(current)

            cand = iso_cand | aniso_cand
            if not np.any(cand):
                break

            current = current | cand

            cur_vals_now = astro_norm_clean[current]
            new_med = float(np.median(cur_vals_now))
            if new_med < med_current * 0.90:
                break
            med_current = new_med
            local_floor = max(seed_p20, med_current * 0.7)

        m = current
        m = binary_fill_holes(m)
        m = binary_closing(m, structure=np.ones((3, 3, 3), dtype=bool))
        m = binary_fill_holes(m)
        m = remove_small_objects(m, min_size=10)

        return m

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

    def measure_cells_2d(self, labels, intensity, nuclei, filename,
                         granules=None, cell_channel=None):
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
        granule_2d = None
        if granules is not None and np.any(granules > 0):
            granule_2d = granules[0] if granules.ndim == 3 else granules

        props = measure.regionprops(labels_2d, intensity_image=intensity_2d)

        if not self.ctx.scale_initialized:
            self.ctx.xy_scale_um = 0.5
            self.ctx.z_scale_um = 2.0

        xy_um = self.ctx.xy_scale_um * self.ctx.downsample
        pixel_area_um2 = xy_um * xy_um

        def _want(name):
            return metric_enabled(self.ctx.metrics_config, name, self.ctx.is_2d)

        need_intensity_stats = _want("cell_intensity_stats") or _want("cell_intensity_percentiles") or \
                               _want("cell_intensity_entropy_gini") or _want("cell_intensity_highlow")
        need_cell_glcm = _want("cell_glcm")
        need_process = _want("cell_process")
        need_snr_cnr_sbr = _want("cell_snr_cnr_sbr")

        need_soma_stats = _want("soma_intensity_stats")
        need_soma_entropy = _want("soma_intensity_entropy")
        need_soma_glcm = _want("soma_glcm")

        need_nuc_stats = _want("nucleus_intensity_stats")
        need_nuc_entropy = _want("nucleus_intensity_entropy")
        need_nuc_glcm = _want("nucleus_glcm")

        data = []
        total = len(props)
        logger = self._make_progress_logger(total, prefix='    ')

        for p in props:
            cid = p.label
            cell_mask_2d = (labels_2d == cid)
            if not np.any(cell_mask_2d):
                logger.update(1)
                continue

            cy, cx = p.centroid
            area_px = float(p.area)
            perim_px = float(getattr(p, 'perimeter', 0.0))
            area_um2 = area_px * pixel_area_um2
            perim_um = perim_px * xy_um

            form_factor = (4.0 * math.pi * area_px / (perim_px ** 2 + 1e-12)) if perim_px > 0 else 0.0
            form_factor = float(max(0.0, min(1.0, form_factor)))

            major = float(getattr(p, 'major_axis_length', 0.0))
            minor = float(getattr(p, 'minor_axis_length', 0.0))
            aspect = (major / minor) if minor > 1e-9 else 1.0
            roundness = (4.0 * area_px / (math.pi * major ** 2 + 1e-12)) if major > 0 else 0.0
            roundness = float(max(0.0, min(1.0, roundness)))
            eccentricity = float(getattr(p, 'eccentricity', 0.0))
            orientation = float(getattr(p, 'orientation', 0.0))

            convex_area = float(getattr(p, 'convex_area', area_px))
            convexity = area_px / (convex_area + 1e-12)
            extent = float(getattr(p, 'extent', 0.0))
            equiv_diam_um = float(getattr(p, 'equivalent_diameter', 0.0)) * xy_um

            minr, minc, maxr, maxc = p.bbox
            bbox_area_um2 = (maxr - minr) * (maxc - minc) * pixel_area_um2
            convex_hull_area_um2 = convex_area * pixel_area_um2

            feret_max_um = float(np.hypot(maxr - minr, maxc - minc)) * xy_um
            feret_min_um = float(min((maxr - minr), (maxc - minc))) * xy_um
            feret_ratio = feret_max_um / (feret_min_um + 1e-12)

            if hasattr(p, 'inertia_tensor_eigenvalues'):
                ev = p.inertia_tensor_eigenvalues
                if len(ev) >= 2 and ev[0] > 0:
                    elongation = float(ev[-1] / ev[0])
                else:
                    elongation = 1.0
            else:
                elongation = 1.0

            cell_vals = intensity_2d[cell_mask_2d]
            cell_mean_int = float(np.mean(cell_vals)) if cell_vals.size > 0 else 0.0

            het = None
            if need_intensity_stats:
                het = self._intensity_heterogeneity_metrics(cell_vals)

            local_var = 0.0
            if _want("cell_local_variance"):
                local_var = self._local_variance_mean(intensity_2d, cell_mask_2d, win=5)

            tex = {'texture_contrast': 0.0, 'texture_homogeneity': 0.0, 'texture_energy': 0.0, 'texture_correlation': 0.0}
            if need_cell_glcm:
                tex = self._glcm_texture(intensity_2d, cell_mask_2d, levels=16)

            fd = 0.0
            if _want("cell_fractal"):
                fd = self._fractal_dimension_box(cell_mask_2d)

            snr = 0.0
            cnr = 0.0
            sbr = 0.0
            if need_snr_cnr_sbr:
                ring_mask = binary_dilation(cell_mask_2d, iterations=3) & (~cell_mask_2d)
                bg_vals = intensity_2d[ring_mask] if np.any(ring_mask) else np.array([], dtype=np.float32)
                bg_mean = float(np.mean(bg_vals)) if bg_vals.size > 0 else 0.0
                bg_std = float(np.std(bg_vals)) if bg_vals.size > 0 else 0.0

                snr = (cell_mean_int - bg_mean) / (bg_std + 1e-12) if bg_vals.size > 0 else 0.0
                cnr = ((cell_mean_int - bg_mean) / (cell_mean_int + bg_mean + 1e-12)) if bg_vals.size > 0 else 0.0
                sbr = cell_mean_int / (bg_mean + 1e-12) if bg_vals.size > 0 else 0.0

            soma_mask_2d = np.zeros_like(cell_mask_2d, dtype=bool)
            if self.ctx._soma_labels is not None:
                if self.ctx._soma_labels.ndim == 3:
                    soma_mask_2d = (self.ctx._soma_labels[0] == cid)
                else:
                    soma_mask_2d = (self.ctx._soma_labels == cid)

            soma_area_px = float(np.sum(soma_mask_2d))
            soma_props = measure.regionprops(soma_mask_2d.astype(np.uint8), intensity_image=intensity_2d) if soma_area_px > 0 else []

            soma_area_um2 = soma_area_px * pixel_area_um2
            soma_perim_um = float(soma_props[0].perimeter) * xy_um if soma_props else 0.0
            soma_solidity = float(soma_props[0].solidity) if soma_props else 1.0
            soma_mean_int = 0.0
            sy, sx = cy, cx

            if soma_area_px > 0:
                soma_int = intensity_2d[soma_mask_2d]
                soma_mean_int = float(np.mean(soma_int))
                sy, sx = soma_props[0].centroid if soma_props else (cy, cx)

            soma_ext = {}
            if soma_area_px > 0 and soma_props:
                sp = soma_props[0]
                if _want("soma_elongation"):
                    soma_ext['soma_elongation'] = float(sp.eccentricity)
                if _want("soma_form_factor"):
                    soma_ext['soma_form_factor'] = float(4.0 * math.pi * sp.area / (sp.perimeter**2 + 1e-12)) if sp.perimeter > 0 else 0.0
                if _want("soma_roundness"):
                    soma_ext['soma_roundness'] = float(4.0 * sp.area / (math.pi * sp.major_axis_length**2 + 1e-12)) if sp.major_axis_length > 0 else 0.0
                if _want("soma_eccentricity"):
                    soma_ext['soma_eccentricity'] = float(sp.eccentricity)
                if _want("soma_convexity"):
                    soma_ext['soma_convexity'] = float(sp.area / (sp.convex_area + 1e-12))

                if need_soma_stats:
                    het_soma = self._intensity_heterogeneity_metrics(soma_int)
                    soma_ext['soma_intensity_std'] = het_soma['intensity_std']
                    soma_ext['soma_intensity_cv'] = het_soma['intensity_cv']
                if need_soma_entropy:
                    het_soma = self._intensity_heterogeneity_metrics(soma_int)
                    soma_ext['soma_intensity_entropy'] = het_soma['intensity_entropy']
                if need_soma_glcm:
                    tex_soma = self._glcm_texture(intensity_2d, soma_mask_2d)
                    soma_ext.update({f"soma_{k}": v for k, v in tex_soma.items()})

            nuc_area_px = 0.0
            nuc_area_um2 = 0.0
            nuc_perim_um = 0.0
            nuc_solidity = 0.0
            nuc_mean_int = 0.0
            ny, nx = cy, cx
            nuc_ext = {}
            nuc_props = []

            if np.any(nuclei_2d > 0):
                nuclei_in_cell = np.unique(nuclei_2d[cell_mask_2d])
                nuclei_in_cell = nuclei_in_cell[nuclei_in_cell > 0]
                if len(nuclei_in_cell) > 0:
                    nuc_mask_2d = (nuclei_2d == nuclei_in_cell[0])
                    nuc_area_px = float(np.sum(nuc_mask_2d))
                    nuc_props = measure.regionprops(nuc_mask_2d.astype(np.uint8), intensity_image=intensity_2d)
                    nuc_area_um2 = nuc_area_px * pixel_area_um2
                    if nuc_props:
                        nuc_int = intensity_2d[nuc_mask_2d]
                        nuc_mean_int = float(np.mean(nuc_int)) if nuc_int.size > 0 else 0.0
                        nuc_perim_um = float(nuc_props[0].perimeter) * xy_um
                        nuc_solidity = float(nuc_props[0].solidity)
                        ny, nx = nuc_props[0].centroid

                        np_ = nuc_props[0]
                        if _want("nucleus_elongation"):
                            nuc_ext['nucleus_elongation'] = float(np_.eccentricity)
                        if _want("nucleus_form_factor"):
                            nuc_ext['nucleus_form_factor'] = float(4.0 * math.pi * np_.area / (np_.perimeter**2 + 1e-12)) if np_.perimeter > 0 else 0.0
                        if _want("nucleus_roundness"):
                            nuc_ext['nucleus_roundness'] = float(4.0 * np_.area / (math.pi * np_.major_axis_length**2 + 1e-12)) if np_.major_axis_length > 0 else 0.0
                        if _want("nucleus_eccentricity"):
                            nuc_ext['nucleus_eccentricity'] = float(np_.eccentricity)
                        if _want("nucleus_convexity"):
                            nuc_ext['nucleus_convexity'] = float(np_.area / (np_.convex_area + 1e-12))

                        if need_nuc_stats:
                            het_nuc = self._intensity_heterogeneity_metrics(nuc_int)
                            nuc_ext['nucleus_intensity_std'] = het_nuc['intensity_std']
                            nuc_ext['nucleus_intensity_cv'] = het_nuc['intensity_cv']
                        if need_nuc_entropy:
                            het_nuc = self._intensity_heterogeneity_metrics(nuc_int)
                            nuc_ext['nucleus_intensity_entropy'] = het_nuc['intensity_entropy']
                        if need_nuc_glcm:
                            tex_nuc = self._glcm_texture(intensity_2d, nuc_mask_2d)
                            nuc_ext.update({f"nucleus_{k}": v for k, v in tex_nuc.items()})

            cyto_mask_2d = cell_mask_2d & (~soma_mask_2d)
            if np.any(cyto_mask_2d):
                cyto_vals = intensity_2d[cyto_mask_2d]
                cyto_mean_int = float(np.mean(cyto_vals)) if cyto_vals.size > 0 else 0.0
            else:
                cyto_mean_int = 0.0
            cyto_area_um2 = max(0.0, area_um2 - soma_area_um2)

            nuc_cyto_ratio = 0.0
            if cyto_area_um2 > 0:
                nuc_cyto_ratio = nuc_area_um2 / cyto_area_um2

            nuc_int_ratio = (nuc_mean_int / (cyto_mean_int + 1e-12)) if cyto_mean_int > 0 else 0.0

            soma_dist = float(np.hypot(sy - cy, sx - cx)) * xy_um
            nuc_dist = float(np.hypot(ny - cy, nx - cx)) * xy_um

            process_count = 0
            process_max_radius = 0
            if need_process and soma_area_px > 0:
                res = self.count_processes_ellipsoid_method(labels, self.ctx._soma_labels)
                process_count = res['process_max_count']
                process_max_radius = res['process_max_radius']

            def _moments(mask):
                ys, xs = np.nonzero(mask)
                if ys.size == 0:
                    return 0.0, 0.0, 0.0, 0.0
                ys = ys.astype(np.float64); xs = xs.astype(np.float64)
                my = float(np.mean(ys)); mx = float(np.mean(xs))
                dy = ys - my; dx = xs - mx
                r2 = dy * dy + dx * dx
                r = np.sqrt(r2)
                if r.sum() <= 0:
                    return 0.0, 0.0, 0.0, 0.0
                w = r / r.sum()
                m1 = float(np.sum(w * r))
                m2 = float(np.sum(w * (r - m1) ** 2))
                m3 = float(np.sum(w * (r - m1) ** 3)) / (m2 ** 1.5 + 1e-12)
                m4 = float(np.sum(w * (r - m1) ** 4)) / (m2 ** 2 + 1e-12)
                return m1, m2, m3, m4

            mom1, mom2, mom3, mom4 = _moments(cell_mask_2d)

            granule_fields = {
                'granule_mean_intensity': 0.0,
                'granule_max_intensity': 0.0,
                'granule_p95_intensity': 0.0,
                'granule_integrated_intensity': 0.0,
                'granule_area_um2': 0.0,
                'granule_area_fraction': 0.0,
                'granule_count': 0,
                'granule_mean_size_um2': 0.0,
                'granule_size_std_um2': 0.0,
                'granule_max_size_um2': 0.0,
                'granule_mean_eccentricity': 0.0,
                'granule_mean_solidity': 0.0,
                'granule_spatial_dispersion_um': 0.0,
                'granule_nearest_neighbor_mean_um': 0.0,
                'granule_clark_evans_index': 0.0,
                'granule_signal_to_background_ratio': 0.0,
            }
            if granule_2d is not None:
                gran_mask = cell_mask_2d & (granule_2d > 0)
                if np.any(gran_mask):
                    garea_px = float(np.sum(gran_mask))
                    granule_fields['granule_area_um2'] = garea_px * pixel_area_um2
                    granule_fields['granule_area_fraction'] = garea_px / (area_px + 1e-12)

                    if cell_channel is not None:
                        gvals = cell_channel[0][gran_mask] if cell_channel.ndim == 3 else cell_channel[gran_mask]
                        if gvals.size > 0:
                            granule_fields['granule_mean_intensity'] = float(np.mean(gvals))
                            granule_fields['granule_max_intensity'] = float(np.max(gvals))
                            granule_fields['granule_p95_intensity'] = float(np.percentile(gvals, 95))
                            granule_fields['granule_integrated_intensity'] = float(np.sum(gvals))

                    lbl_g, n_g = label(gran_mask, structure=np.ones((3, 3)))
                    props_g = measure.regionprops(lbl_g)
                    if props_g:
                        granule_fields['granule_count'] = int(len(props_g))
                        sizes = np.array([pp.area for pp in props_g], dtype=np.float64) * pixel_area_um2
                        granule_fields['granule_mean_size_um2'] = float(np.mean(sizes))
                        granule_fields['granule_size_std_um2'] = float(np.std(sizes))
                        granule_fields['granule_max_size_um2'] = float(np.max(sizes))
                        granule_fields['granule_mean_eccentricity'] = float(
                            np.mean([pp.eccentricity for pp in props_g]))
                        granule_fields['granule_mean_solidity'] = float(
                            np.mean([pp.solidity for pp in props_g]))
                        cents = np.array([pp.centroid for pp in props_g], dtype=np.float64)
                        if len(cents) >= 2:
                            cent = cents.mean(axis=0)
                            disp = np.sqrt(np.mean(np.sum((cents - cent) ** 2, axis=1))) * xy_um
                            granule_fields['granule_spatial_dispersion_um'] = float(disp)
                            dmat = np.linalg.norm(cents[:, None, :] - cents[None, :, :], axis=2)
                            np.fill_diagonal(dmat, np.inf)
                            nn = dmat.min(axis=1) * xy_um
                            nn_mean = float(np.mean(nn))
                            granule_fields['granule_nearest_neighbor_mean_um'] = nn_mean
                            obs_area_um2 = area_um2
                            if obs_area_um2 > 0 and len(cents) > 0:
                                n_pts = float(len(cents))
                                expected = 0.5 * math.sqrt(obs_area_um2 / n_pts)
                                if expected > 0:
                                    granule_fields['granule_clark_evans_index'] = nn_mean / expected

            row = {
                'cell_id': cid,
                'filename': os.path.basename(filename),
                'centroid_y': cy,
                'centroid_x': cx,
                'cell_area_um2': area_um2,
                'soma_area_um2': soma_area_um2,
                'nucleus_area_um2': nuc_area_um2,
                'cytoplasm_area_um2': cyto_area_um2,
                'cell_perimeter_um': perim_um,
                'soma_perimeter_um': soma_perim_um,
                'nucleus_perimeter_um': nuc_perim_um,
                'cell_solidity': float(p.solidity),
                'soma_solidity': soma_solidity,
                'nucleus_solidity': nuc_solidity,
                'cell_elongation': elongation,
                'elongation': elongation,
                'cell_form_factor': form_factor,
                'cell_circularity': form_factor,
                'cell_roundness': roundness,
                'cell_aspect_ratio': aspect,
                'cell_eccentricity': eccentricity,
                'cell_convexity': convexity,
                'cell_extent': extent,
                'cell_equivalent_diameter_um': equiv_diam_um,
                'cell_major_axis_um': major * xy_um,
                'cell_minor_axis_um': minor * xy_um,
                'cell_max_feret_um': feret_max_um,
                'cell_min_feret_um': feret_min_um,
                'cell_feret_ratio': feret_ratio,
                'cell_orientation_deg': float(np.degrees(orientation)),
                'cell_bbox_area_um2': bbox_area_um2,
                'cell_convex_hull_area_um2': convex_hull_area_um2,
                'cell_mean_intensity': cell_mean_int,
                'soma_mean_intensity': soma_mean_int,
                'nucleus_mean_intensity': nuc_mean_int,
                'cytoplasm_mean_intensity': cyto_mean_int,
                'cell_local_variance_mean': local_var,
                'cell_fractal_dimension': fd,
                'cell_snr': snr,
                'cell_cn_r': cnr,
                'cell_signal_to_background_ratio': sbr,
                'cell_moment_1_um': mom1 * xy_um,
                'cell_moment_2_um': mom2 * xy_um * xy_um,
                'cell_moment_3': mom3,
                'cell_moment_4': mom4,
                'cell_nucleus_intensity_ratio': nuc_int_ratio,
                'cell_nucleus_centroid_distance_um': nuc_dist,
                'cell_soma_centroid_distance_um': soma_dist,
                'process_max_count': process_count,
                'process_max_radius': process_max_radius,
                'nuclear_cytoplasmic_ratio': nuc_cyto_ratio,
            }

            if het:
                row.update({
                    'cell_intensity_std': het['intensity_std'],
                    'cell_intensity_cv': het['intensity_cv'],
                    'cell_intensity_iqr': het['intensity_iqr'],
                    'cell_intensity_iqr_rel': het['intensity_iqr_rel'],
                    'cell_intensity_skew': het['intensity_skew'],
                    'cell_intensity_kurt': het['intensity_kurt'],
                    'cell_intensity_entropy': het['intensity_entropy'],
                    'cell_intensity_gini': het['intensity_gini'],
                    'cell_intensity_p10': het['intensity_p10'],
                    'cell_intensity_p50': het['intensity_p50'],
                    'cell_intensity_p90': het['intensity_p90'],
                    'cell_high_intensity_fraction': het['high_intensity_fraction'],
                    'cell_low_intensity_fraction': het['low_intensity_fraction'],
                })

            if need_cell_glcm:
                row.update({
                    'cell_texture_contrast': tex['texture_contrast'],
                    'cell_texture_homogeneity': tex['texture_homogeneity'],
                    'cell_texture_energy': tex['texture_energy'],
                    'cell_texture_correlation': tex['texture_correlation'],
                })

            row.update(granule_fields)
            row.update(soma_ext)
            row.update(nuc_ext)
            data.append(row)
            logger.update(1)

        logger.close()
        df = pd.DataFrame(data)
        return self._filter_metrics(df)

    def measure_cells_3d(self, labels, intensity, nuclei, filename,
                         granules=None, cell_channel=None):
        if labels.shape != intensity.shape:
            self.ctx.log(self.ctx.t(
                "log.measurement_warn_shape",
                fn="measure_cells_3d",
                ls=labels.shape,
                is_=intensity.shape,
            ))

        if self.ctx.is_2d:
            return self.measure_cells_2d(labels, intensity, nuclei, filename,
                                         granules=granules,
                                         cell_channel=cell_channel)

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
        spacing = (z_um, xy_um, xy_um)

        def _want(name):
            return metric_enabled(self.ctx.metrics_config, name, self.ctx.is_2d)

        need_cell_surface = _want("cell_surface") or _want("cell_sphericity") or _want("cell_sav")
        need_cell_axes = _want("cell_axes_3d") or _want("cell_elongation_3d") or _want("cell_flatness")
        need_cell_convex = _want("cell_convexity_3d")
        need_intensity_stats = _want("cell_intensity_stats") or _want("cell_intensity_percentiles") or \
                               _want("cell_intensity_entropy_gini") or _want("cell_intensity_highlow")
        need_cell_glcm = _want("cell_glcm")
        need_cell_process = _want("cell_process")
        need_snr_cnr_sbr = _want("cell_snr_cnr_sbr")

        need_soma_stats = _want("soma_intensity_stats")
        need_soma_entropy = _want("soma_intensity_entropy")
        need_soma_glcm = _want("soma_glcm")
        need_nuc_stats = _want("nucleus_intensity_stats")
        need_nuc_entropy = _want("nucleus_intensity_entropy")
        need_nuc_glcm = _want("nucleus_glcm")

        data = []
        total = len(props)
        logger = self._make_progress_logger(total, prefix='    ')

        for p in props:
            cid = p.label
            cell_mask = (labels == cid)
            if not np.any(cell_mask):
                logger.update(1)
                continue

            cz, cy, cx = p.centroid
            vol_vox = float(p.area)
            vol_um3 = vol_vox * voxel_vol_um3

            cell_vals = intensity[cell_mask]
            cell_mean_int = float(np.mean(cell_vals)) if cell_vals.size > 0 else 0.0

            het_cell = None
            if need_intensity_stats:
                het_cell = self._intensity_heterogeneity_metrics(cell_vals)

            surf_um2 = 0.0
            if need_cell_surface:
                try:
                    mask_f32 = np.ascontiguousarray(cell_mask, dtype=np.float32)
                    verts, faces, normals, values = measure.marching_cubes(
                        mask_f32, level=0.5, spacing=spacing
                    )
                    surf_um2 = float(measure.mesh_surface_area(verts, faces))
                except Exception:
                    surf_um2 = 0.0

            sphericity = 0.0
            if surf_um2 > 0:
                sphericity = (math.pi ** (1.0 / 3.0) * (6.0 * vol_um3) ** (2.0 / 3.0) / surf_um2)
                sphericity = float(max(0.0, min(1.0, sphericity)))
            sa_v = surf_um2 / (vol_um3 + 1e-12)

            major_ax = minor_ax = least_ax = 0.0
            elong3d = 1.0
            ecc3d = 0.0
            flatness = 1.0

            if need_cell_axes:
                coords = np.argwhere(cell_mask).astype(np.float64)
                if len(coords) > 3:
                    coords -= coords.mean(axis=0)
                    coords_phys = coords * np.array([z_um, xy_um, xy_um])
                    inertia = np.dot(coords_phys.T, coords_phys)
                    try:
                        ev = np.linalg.eigvalsh(inertia)
                        ev = np.sort(ev)
                        if ev[2] > 1e-6:
                            major_ax = float(math.sqrt(ev[2]) * 4.0)
                            minor_ax = float(math.sqrt(ev[1]) * 4.0)
                            least_ax = float(math.sqrt(ev[0]) * 4.0)
                            elong3d = float(ev[2] / ev[0]) if ev[0] > 1e-6 else 1.0
                            ecc3d = float(math.sqrt(max(0.0, 1.0 - ev[0] / ev[2])))
                            flatness = float(ev[0] / ev[1]) if ev[1] > 1e-6 else 1.0
                    except Exception:
                        pass

            conv_vol = vol_um3
            if need_cell_convex:
                try:
                    from skimage.morphology import convex_hull_image
                    conv = convex_hull_image(np.ascontiguousarray(cell_mask))
                    conv_vol = float(np.sum(conv)) * voxel_vol_um3
                except Exception:
                    pass

            convexity_3d = vol_um3 / (conv_vol + 1e-12)

            cell_solidity = float(p.solidity) if hasattr(p, 'solidity') else 1.0
            if cell_solidity < 0.05:
                cell_solidity = convexity_3d

            minr, minc, minz, maxr, maxc, maxz = p.bbox
            bbox_vol = max(1, (maxr - minr)) * max(1, (maxc - minc)) * max(1, (maxz - minz)) * voxel_vol_um3
            extent_3d = vol_um3 / (bbox_vol + 1e-12)

            tex_cell = {'texture_contrast': 0.0, 'texture_homogeneity': 0.0, 'texture_energy': 0.0, 'texture_correlation': 0.0}
            if need_cell_glcm:
                slice_mid = int(round(cz))
                slice_mid = max(0, min(labels.shape[0] - 1, slice_mid))
                tex_cell = self._glcm_texture(intensity[slice_mid], cell_mask[slice_mid])

            snr = 0.0
            cnr = 0.0
            sbr = 0.0
            if need_snr_cnr_sbr:
                ring = binary_dilation(cell_mask, iterations=2) & (~cell_mask)
                bg_vals = intensity[ring] if np.any(ring) else np.array([], dtype=np.float32)
                bg_mean = float(np.mean(bg_vals)) if bg_vals.size > 0 else 0.0
                bg_std = float(np.std(bg_vals)) if bg_vals.size > 0 else 0.0

                snr = (cell_mean_int - bg_mean) / (bg_std + 1e-12) if bg_vals.size > 0 else 0.0
                cnr = ((cell_mean_int - bg_mean) / (cell_mean_int + bg_mean + 1e-12)) if bg_vals.size > 0 else 0.0
                sbr = cell_mean_int / (bg_mean + 1e-12) if bg_vals.size > 0 else 0.0

            process_count = 0
            process_max_radius = 0

            soma_mask = (self.ctx._soma_labels == cid)
            soma_vol_vox = float(np.sum(soma_mask))
            soma_props = measure.regionprops(soma_mask.astype(np.uint8), intensity_image=intensity) if soma_vol_vox > 0 else []

            if need_cell_process and soma_vol_vox > 0:
                res = self.count_processes_ellipsoid_method(cell_mask, soma_mask)
                process_count = res['process_max_count']
                process_max_radius = res['process_max_radius']

            soma_ext = {}
            soma_vol_um3 = soma_vol_vox * voxel_vol_um3
            soma_solidity = float(soma_props[0].solidity) if soma_props else 1.0
            soma_mean_int = 0.0
            sy, sx, sz = cy, cx, cz
            soma_sz = cz

            if soma_vol_vox > 0:
                soma_vals = intensity[soma_mask]
                soma_mean_int = float(np.mean(soma_vals))
                sy, sx, sz = (soma_props[0].centroid[::-1] if soma_props else (cy, cx, cz))
                soma_sz = float(soma_props[0].centroid[0]) if soma_props else cz

                if need_soma_stats:
                    het_soma = self._intensity_heterogeneity_metrics(soma_vals)
                    soma_ext['soma_intensity_std'] = het_soma['intensity_std']
                    soma_ext['soma_intensity_cv'] = het_soma['intensity_cv']
                if need_soma_entropy:
                    het_soma = self._intensity_heterogeneity_metrics(soma_vals)
                    soma_ext['soma_intensity_entropy'] = het_soma['intensity_entropy']
                if need_soma_glcm and soma_props:
                    sz_mid = int(round(soma_props[0].centroid[0]))
                    tex_soma = self._glcm_texture(intensity[sz_mid], soma_mask[sz_mid])
                    soma_ext.update({f"soma_{k}": v for k, v in tex_soma.items()})

            nuc_vol_vox = 0.0
            nuc_vol_um3 = 0.0
            nuc_solidity = 0.0
            nuc_mean_int = 0.0
            nz_c, ny_c, nx_c = cz, cy, cx
            nuc_ext = {}

            if np.any(nuclei > 0):
                nuclei_in_cell = np.unique(nuclei[cell_mask])
                nuclei_in_cell = nuclei_in_cell[nuclei_in_cell > 0]
                if len(nuclei_in_cell) > 0:
                    nuc_mask = (nuclei == nuclei_in_cell[0])
                    nuc_vol_vox = float(np.sum(nuc_mask))
                    nuc_props = measure.regionprops(nuc_mask.astype(np.uint8), intensity_image=intensity)
                    nuc_vol_um3 = nuc_vol_vox * voxel_vol_um3

                    if nuc_props:
                        nuc_vals = intensity[nuc_mask]
                        nuc_mean_int = float(np.mean(nuc_vals)) if nuc_vals.size > 0 else 0.0
                        nuc_solidity = float(nuc_props[0].solidity)
                        nz_c, ny_c, nx_c = nuc_props[0].centroid

                        if need_nuc_stats:
                            het_nuc = self._intensity_heterogeneity_metrics(nuc_vals)
                            nuc_ext['nucleus_intensity_std'] = het_nuc['intensity_std']
                            nuc_ext['nucleus_intensity_cv'] = het_nuc['intensity_cv']
                        if need_nuc_entropy:
                            het_nuc = self._intensity_heterogeneity_metrics(nuc_vals)
                            nuc_ext['nucleus_intensity_entropy'] = het_nuc['intensity_entropy']
                        if need_nuc_glcm:
                            nz_mid = int(round(nuc_props[0].centroid[0]))
                            tex_nuc = self._glcm_texture(intensity[nz_mid], nuc_mask[nz_mid])
                            nuc_ext.update({f"nucleus_{k}": v for k, v in tex_nuc.items()})

            cyto_mask = cell_mask & (~soma_mask)
            if np.any(cyto_mask):
                cyto_vals = intensity[cyto_mask]
                cyto_mean_int = float(np.mean(cyto_vals)) if cyto_vals.size > 0 else 0.0
            else:
                cyto_mean_int = 0.0

            cyto_vol_um3 = max(0.0, vol_um3 - nuc_vol_um3)
            nuc_cyto_ratio = (nuc_vol_um3 / cyto_vol_um3) if cyto_vol_um3 > 0 else 0.0
            nuc_int_ratio = (nuc_mean_int / (cyto_mean_int + 1e-12)) if cyto_mean_int > 0 else 0.0

            soma_dist = float(np.sqrt(((sy - cy) * xy_um) ** 2 + ((sx - cx) * xy_um) ** 2 + ((soma_sz - cz) * z_um) ** 2))
            nuc_dist = float(np.sqrt(((ny_c - cy) * xy_um) ** 2 + ((nx_c - cx) * xy_um) ** 2 + ((nz_c - cz) * z_um) ** 2))

            row = {
                'cell_id': cid,
                'filename': os.path.basename(filename),
                'centroid_z': cz,
                'centroid_y': cy,
                'centroid_x': cx,
                'cell_volume_um3': vol_um3,
                'soma_volume_um3': soma_vol_um3,
                'nucleus_volume_um3': nuc_vol_um3,
                'cytoplasm_volume_um3': cyto_vol_um3,
                'cell_surface_area_um2': surf_um2,
                'cell_surface_to_volume_ratio': sa_v,
                'cell_sphericity': sphericity,
                'cell_major_axis_um': major_ax,
                'cell_minor_axis_um': minor_ax,
                'cell_least_axis_um': least_ax,
                'cell_elongation_3d': elong3d,
                'cell_eccentricity': ecc3d,
                'cell_flatness': flatness,
                'cell_convex_volume_um3': conv_vol,
                'cell_convexity_3d': convexity_3d,
                'cell_extent_3d': extent_3d,
                'cell_solidity': cell_solidity,
                'soma_solidity': soma_solidity,
                'nucleus_solidity': nuc_solidity,
                'elongation': elong3d,
                'cell_mean_intensity': cell_mean_int,
                'soma_mean_intensity': soma_mean_int,
                'nucleus_mean_intensity': nuc_mean_int,
                'cytoplasm_mean_intensity': cyto_mean_int,
                'cell_signal_to_background_ratio': sbr,
                'cell_snr': snr,
                'cell_cn_r': cnr,
                'cell_nucleus_intensity_ratio': nuc_int_ratio,
                'cell_nucleus_centroid_distance_um': nuc_dist,
                'cell_soma_centroid_distance_um': soma_dist,
                'process_max_count': process_count,
                'process_max_radius': process_max_radius,
                'nuclear_cytoplasmic_ratio': nuc_cyto_ratio,
            }

            if het_cell:
                row.update({
                    'cell_intensity_std': het_cell['intensity_std'],
                    'cell_intensity_cv': het_cell['intensity_cv'],
                    'cell_intensity_iqr': het_cell['intensity_iqr'],
                    'cell_intensity_iqr_rel': het_cell['intensity_iqr_rel'],
                    'cell_intensity_skew': het_cell['intensity_skew'],
                    'cell_intensity_kurt': het_cell['intensity_kurt'],
                    'cell_intensity_entropy': het_cell['intensity_entropy'],
                    'cell_intensity_gini': het_cell['intensity_gini'],
                    'cell_intensity_p10': het_cell['intensity_p10'],
                    'cell_intensity_p50': het_cell['intensity_p50'],
                    'cell_intensity_p90': het_cell['intensity_p90'],
                    'cell_high_intensity_fraction': het_cell['high_intensity_fraction'],
                    'cell_low_intensity_fraction': het_cell['low_intensity_fraction'],
                })

            if need_cell_glcm:
                row.update({
                    'cell_texture_contrast': tex_cell['texture_contrast'],
                    'cell_texture_homogeneity': tex_cell['texture_homogeneity'],
                    'cell_texture_energy': tex_cell['texture_energy'],
                    'cell_texture_correlation': tex_cell['texture_correlation'],
                })

            row.update(soma_ext)
            row.update(nuc_ext)
            data.append(row)
            logger.update(1)

        logger.close()
        df = pd.DataFrame(data)
        return self._filter_metrics(df)

    def _intensity_heterogeneity_metrics(self, values):
        v = np.asarray(values, dtype=np.float32).ravel()
        v = v[np.isfinite(v)]
        if v.size < 2:
            return {
                'intensity_std': 0.0, 'intensity_cv': 0.0,
                'intensity_iqr': 0.0, 'intensity_iqr_rel': 0.0,
                'intensity_skew': 0.0, 'intensity_kurt': 0.0,
                'intensity_entropy': 0.0, 'intensity_gini': 0.0,
                'intensity_p10': 0.0, 'intensity_p50': 0.0, 'intensity_p90': 0.0,
                'high_intensity_fraction': 0.0, 'low_intensity_fraction': 0.0,
            }
        mean_v = float(np.mean(v))
        std_v = float(np.std(v))
        med_v = float(np.median(v))
        p10, p25, p50, p75, p90 = np.percentile(v, [10, 25, 50, 75, 90])
        iqr = float(p75 - p25)

        m3 = float(np.mean((v - mean_v) ** 3))
        m4 = float(np.mean((v - mean_v) ** 4))
        skew = m3 / (std_v ** 3 + 1e-12) if std_v > 0 else 0.0
        kurt = m4 / (std_v ** 4 + 1e-12) - 3.0 if std_v > 0 else 0.0

        hist, _ = np.histogram(v, bins=32, density=False)
        hist = hist.astype(np.float64)
        hist_sum = hist.sum()
        if hist_sum > 0:
            p = hist / hist_sum
            p = p[p > 0]
            entropy = float(-np.sum(p * np.log2(p)))
        else:
            entropy = 0.0

        v_sorted = np.sort(v.astype(np.float64))
        n = v_sorted.size
        cum = np.cumsum(v_sorted)
        total = cum[-1] if cum[-1] > 0 else 1e-12
        gini = float((n + 1 - 2 * np.sum(cum) / total) / n)
        gini = max(0.0, min(1.0, gini))

        high_frac = float(np.mean(v > p75))
        low_frac = float(np.mean(v < p25))

        return {
            'intensity_std': std_v,
            'intensity_cv': std_v / (mean_v + 1e-12),
            'intensity_iqr': iqr,
            'intensity_iqr_rel': iqr / (med_v + 1e-12),
            'intensity_skew': float(skew),
            'intensity_kurt': float(kurt),
            'intensity_entropy': entropy,
            'intensity_gini': gini,
            'intensity_p10': float(p10),
            'intensity_p50': float(p50),
            'intensity_p90': float(p90),
            'high_intensity_fraction': high_frac,
            'low_intensity_fraction': low_frac,
        }

    def _local_variance_mean(self, image, mask, win=5):
        if not np.any(mask):
            return 0.0
        img = image.astype(np.float32)
        m = mask.astype(np.float32)
        k = np.ones((win, win), dtype=np.float32)
        k /= k.sum()
        mean_local = convolve(img * m, k, mode='constant')
        mean_sq = convolve((img ** 2) * m, k, mode='constant')
        var_local = np.clip(mean_sq - mean_local ** 2, 0, None)
        vals = var_local[mask > 0]
        return float(np.mean(vals)) if vals.size > 0 else 0.0

    def _glcm_texture(self, image, mask, levels=16):
        if not np.any(mask):
            return {'texture_contrast': 0.0, 'texture_homogeneity': 0.0,
                    'texture_energy': 0.0, 'texture_correlation': 0.0}
        img = image.astype(np.float32)
        vals = img[mask]
        if vals.size < 4:
            return {'texture_contrast': 0.0, 'texture_homogeneity': 0.0,
                    'texture_energy': 0.0, 'texture_correlation': 0.0}
        lo, hi = float(np.min(vals)), float(np.max(vals))
        if hi - lo < 1e-9:
            return {'texture_contrast': 0.0, 'texture_homogeneity': 0.0,
                    'texture_energy': 0.0, 'texture_correlation': 0.0}
        q = np.clip(((img - lo) / (hi - lo) * (levels - 1)).astype(np.int32), 0, levels - 1)

        directions = [(0, 1), (1, 0), (1, 1), (1, -1)]
        h, w = q.shape
        acc = np.zeros((levels, levels), dtype=np.float64)
        for dy, dx in directions:
            ys0 = max(0, -dy); ys1 = h - max(0, dy)
            xs0 = max(0, -dx); xs1 = w - max(0, dx)
            if ys1 <= ys0 or xs1 <= xs0:
                continue
            a = q[ys0:ys1, xs0:xs1]
            b = q[ys0 + dy:ys1 + dy, xs0 + dx:xs1 + dx]
            ma = mask[ys0:ys1, xs0:xs1]
            mb = mask[ys0 + dy:ys1 + dy, xs0 + dx:xs1 + dx]
            valid = ma & mb
            if not np.any(valid):
                continue
            aa = a[valid].ravel(); bb = b[valid].ravel()
            np.add.at(acc, (aa, bb), 1.0)
            np.add.at(acc, (bb, aa), 1.0)
        s = acc.sum()
        if s <= 0:
            return {'texture_contrast': 0.0, 'texture_homogeneity': 0.0,
                    'texture_energy': 0.0, 'texture_correlation': 0.0}
        p = acc / s
        i_idx = np.arange(levels).reshape(-1, 1)
        j_idx = np.arange(levels).reshape(1, -1)
        diff = i_idx - j_idx
        contrast = float(np.sum((diff ** 2) * p))
        homogeneity = float(np.sum(p / (1.0 + diff ** 2)))
        energy = float(np.sqrt(np.sum(p ** 2)))
        mu_i = float(np.sum(i_idx * p))
        mu_j = float(np.sum(j_idx * p))
        var_i = float(np.sum((i_idx - mu_i) ** 2 * p))
        var_j = float(np.sum((j_idx - mu_j) ** 2 * p))
        if var_i > 0 and var_j > 0:
            corr = float(np.sum((i_idx - mu_i) * (j_idx - mu_j) * p) / np.sqrt(var_i * var_j))
        else:
            corr = 0.0
        return {'texture_contrast': contrast, 'texture_homogeneity': homogeneity,
                'texture_energy': energy, 'texture_correlation': corr}

    def _fractal_dimension_box(self, mask_2d):
        m = mask_2d.astype(bool)
        if m.sum() < 10:
            return 0.0
        sizes = [2, 4, 8, 16, 32]
        counts = []
        h, w = m.shape
        for s in sizes:
            ys = np.arange(0, h, s)
            xs = np.arange(0, w, s)
            n = 0
            for y in ys:
                for x in xs:
                    if np.any(m[y:y + s, x:x + s]):
                        n += 1
            counts.append(max(n, 1))
        log_s = np.log(1.0 / np.array(sizes, dtype=np.float64))
        log_c = np.log(np.array(counts, dtype=np.float64))
        A = np.vstack([log_s, np.ones_like(log_s)]).T
        try:
            slope, _ = np.linalg.lstsq(A, log_c, rcond=None)[0]
            return float(slope)
        except Exception:
            return 0.0

    def _filter_metrics(self, df):
        cfg = getattr(self.ctx, "metrics_config", None)
        if not cfg or df is None or len(df) == 0:
            return df

        cell_cols = {
            "cell_area_um2": "cell_area",
            "cell_perimeter_um": "cell_perimeter",
            "cell_solidity": "cell_solidity",
            "cell_elongation": "cell_elongation",
            "elongation": "cell_elongation",
            "cell_form_factor": "cell_form_factor",
            "cell_circularity": "cell_form_factor",
            "cell_roundness": "cell_roundness",
            "cell_aspect_ratio": "cell_aspect_ratio",
            "cell_eccentricity": "cell_eccentricity",
            "cell_convexity": "cell_convexity",
            "cell_extent": "cell_extent",
            "cell_major_axis_um": "cell_major_minor",
            "cell_minor_axis_um": "cell_major_minor",
            "cell_max_feret_um": "cell_feret",
            "cell_min_feret_um": "cell_feret",
            "cell_feret_ratio": "cell_feret",
            "cell_bbox_area_um2": "cell_bbox_hull",
            "cell_convex_hull_area_um2": "cell_bbox_hull",
            "cell_orientation_deg": "cell_bbox_hull",
            "cell_equivalent_diameter_um": "cell_equiv_diam",
            "cell_moment_1_um": "cell_moments",
            "cell_moment_2_um": "cell_moments",
            "cell_moment_3": "cell_moments",
            "cell_moment_4": "cell_moments",
            "cell_volume_um3": "cell_volume",
            "cell_surface_area_um2": "cell_surface",
            "cell_sphericity": "cell_sphericity",
            "cell_surface_to_volume_ratio": "cell_sav",
            "cell_convex_volume_um3": "cell_convexity_3d",
            "cell_convexity_3d": "cell_convexity_3d",
            "cell_extent_3d": "cell_extent_3d",
            "cell_flatness": "cell_axes_3d",
            "cell_least_axis_um": "cell_axes_3d",
            "cell_elongation_3d": "cell_axes_3d",
            "cell_mean_intensity": "cell_mean_intensity",
            "cytoplasm_mean_intensity": "cell_mean_intensity",
            "cell_intensity_std": "cell_intensity_stats",
            "cell_intensity_cv": "cell_intensity_stats",
            "cell_intensity_iqr": "cell_intensity_stats",
            "cell_intensity_iqr_rel": "cell_intensity_stats",
            "cell_intensity_skew": "cell_intensity_stats",
            "cell_intensity_kurt": "cell_intensity_stats",
            "cell_intensity_p10": "cell_intensity_percentiles",
            "cell_intensity_p50": "cell_intensity_percentiles",
            "cell_intensity_p90": "cell_intensity_percentiles",
            "cell_intensity_entropy": "cell_intensity_entropy_gini",
            "cell_intensity_gini": "cell_intensity_entropy_gini",
            "cell_high_intensity_fraction": "cell_intensity_highlow",
            "cell_low_intensity_fraction": "cell_intensity_highlow",
            "cell_snr": "cell_snr_cnr_sbr",
            "cell_cn_r": "cell_snr_cnr_sbr",
            "cell_signal_to_background_ratio": "cell_snr_cnr_sbr",
            "cell_local_variance_mean": "cell_local_variance",
            "cell_texture_contrast": "cell_glcm",
            "cell_texture_homogeneity": "cell_glcm",
            "cell_texture_energy": "cell_glcm",
            "cell_texture_correlation": "cell_glcm",
            "cell_fractal_dimension": "cell_fractal",
            "process_max_count": "cell_process",
            "process_max_radius": "cell_process",
            "cell_skeleton_length_um": "cell_skeleton",
            "cell_skeleton_branch_points": "cell_skeleton",
            "cell_skeleton_endpoints": "cell_skeleton",
            "cell_skeleton_tortuosity": "cell_skeleton",
        }

        soma_cols = {
            "soma_area_um2": "soma_area",
            "soma_perimeter_um": "soma_perimeter",
            "soma_solidity": "soma_solidity",
            "soma_elongation": "soma_elongation",
            "soma_form_factor": "soma_form_factor",
            "soma_roundness": "soma_roundness",
            "soma_eccentricity": "soma_eccentricity",
            "soma_convexity": "soma_convexity",
            "soma_volume_um3": "soma_volume",
            "soma_elongation_3d": "soma_elongation_3d",
            "soma_flatness": "soma_flatness",
            "soma_extent_3d": "soma_extent_3d",
            "soma_mean_intensity": "soma_mean_intensity",
            "soma_intensity_std": "soma_intensity_stats",
            "soma_intensity_cv": "soma_intensity_stats",
            "soma_intensity_entropy": "soma_intensity_entropy",
            "soma_texture_contrast": "soma_glcm",
            "soma_texture_homogeneity": "soma_glcm",
            "soma_texture_energy": "soma_glcm",
            "soma_texture_correlation": "soma_glcm",
        }

        nucleus_cols = {
            "nucleus_area_um2": "nucleus_area",
            "nucleus_perimeter_um": "nucleus_perimeter",
            "nucleus_solidity": "nucleus_solidity",
            "nucleus_elongation": "nucleus_elongation",
            "nucleus_form_factor": "nucleus_form_factor",
            "nucleus_roundness": "nucleus_roundness",
            "nucleus_eccentricity": "nucleus_eccentricity",
            "nucleus_convexity": "nucleus_convexity",
            "nucleus_volume_um3": "nucleus_volume",
            "nucleus_elongation_3d": "nucleus_elongation_3d",
            "nucleus_flatness": "nucleus_flatness",
            "nucleus_extent_3d": "nucleus_extent_3d",
            "nucleus_mean_intensity": "nucleus_mean_intensity",
            "nucleus_intensity_std": "nucleus_intensity_stats",
            "nucleus_intensity_cv": "nucleus_intensity_stats",
            "nucleus_intensity_entropy": "nucleus_intensity_entropy",
            "nucleus_texture_contrast": "nucleus_glcm",
            "nucleus_texture_homogeneity": "nucleus_glcm",
            "nucleus_texture_energy": "nucleus_glcm",
            "nucleus_texture_correlation": "nucleus_glcm",
        }

        nucleus_group_cols = {
            "cytoplasm_area_um2": "cell_nuc_cyto",
            "cytoplasm_volume_um3": "cell_nuc_cyto",
            "nuclear_cytoplasmic_ratio": "cell_nuc_cyto",
            "cell_nucleus_intensity_ratio": "cell_nuc_intensity",
            "cell_nucleus_centroid_distance_um": "cell_nuc_intensity",
            "cell_soma_centroid_distance_um": "cell_nuc_intensity",
        }

        granule_cols = {
            "granule_mean_intensity": "granule_intensity",
            "granule_max_intensity": "granule_intensity",
            "granule_p95_intensity": "granule_intensity",
            "granule_integrated_intensity": "granule_intensity",
            "granule_area_um2": "granule_size_count",
            "granule_volume_um3": "granule_size_count",
            "granule_area_fraction": "granule_size_count",
            "granule_volume_fraction": "granule_size_count",
            "granule_count": "granule_size_count",
            "granule_mean_size_um2": "granule_size_count",
            "granule_mean_size_um3": "granule_size_count",
            "granule_size_std_um2": "granule_size_count",
            "granule_size_std_um3": "granule_size_count",
            "granule_max_size_um2": "granule_size_count",
            "granule_max_size_um3": "granule_size_count",
            "granule_mean_eccentricity": "granule_shape",
            "granule_mean_solidity": "granule_shape",
            "granule_spatial_dispersion_um": "granule_spatial",
            "granule_nearest_neighbor_mean_um": "granule_spatial",
            "granule_clark_evans_index": "granule_spatial",
            "granule_signal_to_background_ratio": "granule_sbr",
        }

        drop = set()

        for col, key in cell_cols.items():
            if not bool(cfg.get(key, True)) and col in df.columns:
                drop.add(col)

        for col, key in soma_cols.items():
            if not bool(cfg.get(key, True)) and col in df.columns:
                drop.add(col)

        for col, key in nucleus_cols.items():
            if not bool(cfg.get(key, True)) and col in df.columns:
                drop.add(col)

        for col, key in nucleus_group_cols.items():
            if not bool(cfg.get(key, True)) and col in df.columns:
                drop.add(col)

        for col, key in granule_cols.items():
            if not bool(cfg.get(key, True)) and col in df.columns:
                drop.add(col)

        if not getattr(self.ctx, "has_granules", False):
            for c in df.columns:
                if c.startswith("granule_"):
                    drop.add(c)

        if getattr(self.ctx, "is_2d", False):
            for c in df.columns:
                if c in ("cell_volume_um3", "cell_surface_area_um2", "cell_sphericity",
                         "cell_surface_to_volume_ratio", "cell_convex_volume_um3",
                         "cell_convexity_3d", "cell_extent_3d", "cell_flatness",
                         "cell_least_axis_um", "cell_elongation_3d",
                         "soma_volume_um3", "soma_elongation_3d", "soma_flatness",
                         "soma_extent_3d", "soma_convexity_3d",
                         "nucleus_volume_um3", "nucleus_elongation_3d", "nucleus_flatness",
                         "nucleus_extent_3d", "nucleus_convexity_3d",
                         "cytoplasm_volume_um3"):
                    drop.add(c)
        else:
            for c in df.columns:
                if c in ("cell_perimeter_um", "cell_major_axis_um", "cell_minor_axis_um",
                         "cell_max_feret_um", "cell_min_feret_um", "cell_feret_ratio",
                         "cell_bbox_area_um2", "cell_convex_hull_area_um2",
                         "cell_orientation_deg", "cell_equivalent_diameter_um",
                         "cell_moment_1_um", "cell_moment_2_um", "cell_moment_3", "cell_moment_4",
                         "soma_perimeter_um", "soma_elongation", "soma_form_factor",
                         "soma_roundness", "soma_eccentricity", "soma_convexity",
                         "nucleus_perimeter_um", "nucleus_elongation", "nucleus_form_factor",
                         "nucleus_roundness", "nucleus_eccentricity", "nucleus_convexity",
                         "cytoplasm_area_um2"):
                    drop.add(c)

        keep_always = {
            "cell_id", "filename",
            "centroid_y", "centroid_x", "centroid_z",
        }
        drop = drop - keep_always
        return df.drop(columns=[c for c in drop if c in df.columns],
                       errors="ignore")

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
                labeled_intersect, _ = label(intersection,
                                             structure=np.ones((3, 3, 3)))
                props = measure.regionprops(labeled_intersect)
                valid_components = 0
                for prop in props:
                    if prop.area >= 5:
                        if hasattr(prop, 'inertia_tensor_eigenvalues'):
                            ev = prop.inertia_tensor_eigenvalues
                            if len(ev) >= 3 and ev[0] > 0:
                                elongation = ev[2] / ev[0]
                                if elongation > 1.5:
                                    valid_components += 1
                            else:
                                valid_components += 1
                        else:
                            valid_components += 1
                process_counts.append({'count': valid_components})
            else:
                process_counts.append({'count': 0})

        if len(process_counts) > 0:
            counts = [p['count'] for p in process_counts]
            max_processes = max(counts)
            max_idx = np.argmax(counts)
            max_radius = radii_xy[max_idx] if max_idx < len(radii_xy) else 0
            return {
                'process_max_count': max_processes,
                'process_max_radius': max_radius,
            }
        else:
            return {
                'process_max_count': 0,
                'process_max_radius': 0,
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

    def _append_summary_count(self, filepath, cells_count):
        suffix = self._get_cell_suffix()
        dim_suffix = self._get_dim_suffix()
        out_dir = os.path.join(os.path.dirname(filepath), f"all_results_{suffix}_{dim_suffix}")
        os.makedirs(out_dir, exist_ok=True)

        summary_csv = os.path.join(out_dir, "summary_cell_counts.csv")
        current_filename = os.path.basename(filepath)

        new_row = pd.DataFrame([{
            'filename': current_filename,
            'cells_count': int(cells_count)
        }])

        if os.path.exists(summary_csv):
            try:
                existing = pd.read_csv(summary_csv)
            except Exception as e:
                self.ctx.log(self.ctx.t("log.results_read_failed",
                                        path=summary_csv, err=e))
                existing = pd.DataFrame()

            if len(existing) > 0 and 'filename' in existing.columns:
                existing = existing[existing['filename'] != current_filename]

            combined = pd.concat([existing, new_row], ignore_index=True)
        else:
            combined = new_row

        combined.to_csv(summary_csv, index=False)

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

    def save_3d_stacks(self, astro, dapi, base_name, out_dir,
                       soma=None, nuclei=None, cells=None):
        stacks_dir = os.path.join(out_dir, "3d_stacks")
        os.makedirs(stacks_dir, exist_ok=True)

        astro_min = float(astro.min())
        astro_max = float(astro.max())
        astro_norm = ((astro - astro_min) /
                      (astro_max - astro_min + 1e-6) * 65535).astype(np.uint16)
        tifffile.imwrite(
            os.path.join(stacks_dir, f"{base_name}_astro_3d.tiff"),
            astro_norm,
            compression='zlib'
        )

        if dapi is not None and dapi.size > 0:
            dapi_min = float(dapi.min())
            dapi_max = float(dapi.max())
            dapi_norm = ((dapi - dapi_min) /
                         (dapi_max - dapi_min + 1e-6) * 65535).astype(np.uint16)
            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_dapi_3d.tiff"),
                dapi_norm,
                compression='zlib'
            )

        if soma is not None and soma.max() > 0:
            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_soma_ids.tiff"),
                soma.astype(np.uint16),
                compression='zlib'
            )

        if nuclei is not None and nuclei.max() > 0:
            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_nuclei_ids.tiff"),
                nuclei.astype(np.uint16),
                compression='zlib'
            )

        if cells is not None and cells.max() > 0:
            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_cells_ids.tiff"),
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
                rgb_stack[mask] = [int(r * 255), int(g * 255), int(b * 255)]

            tifffile.imwrite(
                os.path.join(stacks_dir, f"{base_name}_cells_rgb.tiff"),
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

    def _save_deconvolved_stacks(self, astro, dapi, base_name, out_dir, granules=None):
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

        if dapi is not None and np.any(dapi > 0):
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

        if granules is not None and granules.size > 0:
            gran_norm = granules.astype(np.float32)
            gmax = np.max(gran_norm)
            if gmax > 0:
                gran_u16 = (gran_norm / gmax * 65535).astype(np.uint16)
            else:
                gran_u16 = gran_norm.astype(np.uint16)

            tifffile.imwrite(
                os.path.join(deconv_dir, f"{base_name}_granules_deconvolved.tiff"),
                gran_u16,
                compression='zlib'
            )
            self.ctx.log(self.ctx.t("log.viz_deconv_granules_saved", shape=granules.shape))

            try:
                mip = np.max(gran_norm, axis=0)
                nzm = mip[mip > 0]
                if len(nzm) > 0:
                    p0, p99 = np.percentile(nzm, [1, 99])
                    mip_norm = np.clip((mip - p0) / (p99 - p0 + 1e-6), 0, 1)
                else:
                    mip_norm = mip

                fig, ax = plt.subplots(1, 1, figsize=(8, 8))
                ax.imshow(mip_norm, cmap='gray')
                ax.set_title(self.ctx.t("log.deconv_granules_title", name=base_name), fontsize=12)
                ax.axis('off')
                plt.tight_layout()
                plt.savefig(
                    os.path.join(deconv_dir, f"{base_name}_granules_deconvolved_mip.png"),
                    dpi=200, bbox_inches='tight'
                )
                plt.close(fig)
            except Exception as e:
                self.ctx.log(self.ctx.t("log.viz_mip_granules_failed", err=e))

class Analyzer:
    def __init__(self):
        self.ctx = AnalysisContext()
        self.file_io = FileIO(self.ctx)
        self.optics = Optics(self.ctx, self.file_io)
        self.deconv = Deconvolution(self.ctx)
        self.segmenter = Segmentation(self.ctx)
        self.measurer = Measurements(self.ctx)
        self.viz = Visualization(self.ctx)
        self.ctx.measurer = self.measurer

    def batch_process_files(self, files):
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

            if not self.ctx._preprocess_decided:
                probe = self.file_io._read_image_universal(filepath)
                probe_y, probe_x = 0, 0
                probe_mip = None
                if probe is not None:
                    if probe.ndim == 2:
                        probe_y, probe_x = probe.shape
                        probe_stack = probe[np.newaxis, ...]
                    elif probe.ndim == 3:
                        if probe.shape[0] <= 8 and probe.shape[1] > probe.shape[0]:
                            probe_y, probe_x = probe.shape[1], probe.shape[2]
                            ch0 = int(self.ctx.cell_channel_idx or 0)
                            ch0 = max(0, min(probe.shape[0] - 1, ch0))
                            probe_stack = probe[ch0:ch0 + 1]
                        else:
                            probe_y, probe_x = probe.shape[1], probe.shape[2]
                            probe_stack = probe[np.newaxis, ...]
                    elif probe.ndim == 4:
                        if probe.shape[0] <= 8:
                            probe_y, probe_x = probe.shape[2], probe.shape[3]
                            ch0 = int(self.ctx.cell_channel_idx or 0)
                            ch0 = max(0, min(probe.shape[0] - 1, ch0))
                            probe_stack = probe[ch0:ch0 + 1]
                        else:
                            probe = np.transpose(probe, (1, 0, 2, 3))
                            probe_y, probe_x = probe.shape[2], probe.shape[3]
                            ch0 = int(self.ctx.cell_channel_idx or 0)
                            ch0 = max(0, min(probe.shape[0] - 1, ch0))
                            probe_stack = probe[ch0:ch0 + 1]
                    elif probe.ndim == 5:
                        probe = probe[0]
                        if probe.shape[0] > 8 and probe.shape[1] <= 8:
                            probe = np.transpose(probe, (1, 0, 2, 3))
                        probe_y, probe_x = probe.shape[2], probe.shape[3]
                        ch0 = int(self.ctx.cell_channel_idx or 0)
                        ch0 = max(0, min(probe.shape[0] - 1, ch0))
                        probe_stack = probe[ch0:ch0 + 1]

                    if probe_stack is not None and probe_stack.size > 0:
                        ps = probe_stack.astype(np.float32)
                        ps_nz = ps[ps > 0]
                        if ps_nz.size > 0:
                            p0, p99 = np.percentile(ps_nz, [1, 99])
                            ps_n = np.clip((ps - p0) / (p99 - p0 + 1e-6), 0, 1)
                        else:
                            ps_n = ps
                        probe_mip = np.max(ps_n, axis=0)
                        while probe_mip.ndim > 2:
                            probe_mip = np.max(probe_mip, axis=0)

                need_ds_q = (max(probe_y, probe_x) > 512)
                xy_um = float(self.ctx.xy_scale_um or 0.0)
                size_um = max(probe_y, probe_x) * xy_um
                need_crop_q = (size_um > 350.0)

                allow_ds = True
                crop_box = None

                if (need_ds_q or need_crop_q) and probe_mip is not None:
                    res = self.ctx.bridge._ask(
                        "preprocess_choice",
                        need_downsample_q=bool(need_ds_q),
                        need_crop_q=bool(need_crop_q),
                        img_w=int(probe_x),
                        img_h=int(probe_y),
                        allow_downsample_default=True,
                        do_crop_default=False,
                    )
                    if res is not None:
                        allow_ds = bool(res.get("allow_downsample", True))
                        do_crop = bool(res.get("do_crop", False))
                        if do_crop:
                            crop_res = self.ctx.bridge._ask(
                                "crop_editor",
                                mip=probe_mip,
                                xy_um=xy_um,
                            )
                            if crop_res is not None:
                                crop_box = tuple(int(v) for v in crop_res)

                if not need_ds_q:
                    allow_ds = True

                self.ctx.allow_downsample = allow_ds
                self.ctx.crop_box = crop_box
                self.ctx._preprocess_decided = True

                if self.ctx.downsample == 1 and allow_ds:
                    self.ctx.downsample = compute_downsample(
                        probe_y, probe_x, target=512
                    ) if need_ds_q else 1
                elif not allow_ds:
                    self.ctx.downsample = 1

                base_probe = os.path.splitext(os.path.basename(filepath))[0]
                suffix_probe = self.measurer._get_cell_suffix()
                pre_dir = os.path.join(
                    os.path.dirname(filepath),
                    f"{base_probe}_{suffix_probe}_analysis"
                )
                os.makedirs(pre_dir, exist_ok=True)
                crop_meta_path = os.path.join(pre_dir, "crop_info.txt")
                try:
                    with open(crop_meta_path, "w", encoding="utf-8") as f:
                        f.write(f"allow_downsample={allow_ds}\n")
                        if crop_box is not None:
                            f.write(
                                f"crop_y0={crop_box[0]}\n"
                                f"crop_y1={crop_box[1]}\n"
                                f"crop_x0={crop_box[2]}\n"
                                f"crop_x1={crop_box[3]}\n"
                            )
                        else:
                            f.write("crop_box=None\n")
                except Exception:
                    pass

            img = self.file_io.load_and_downsample(filepath)
            if img is None:
                return None, None

            astro = None
            dapi = None
            granules_stack = None

            for i, ch_idx in enumerate(self.ctx.selected_channels):
                if ch_idx == self.ctx.cell_channel_idx:
                    astro = img[i]
                elif self.ctx.has_dapi and ch_idx == self.ctx.dapi_channel_idx:
                    dapi = img[i]
                elif self.ctx.has_granules and ch_idx == self.ctx.granule_channel_idx:
                    granules_stack = img[i]

            if astro is None:
                astro = img[0]

            if dapi is None:
                dapi = np.zeros_like(astro)

            if self.ctx.has_dapi and not np.any(dapi > 0):
                self.ctx.log(self.ctx.t("log.dapi_not_used_msg"))

            if self.ctx.deconvolve_images:
                pixel_xy_um = self.ctx.xy_scale_um * self.ctx.downsample
                pixel_z_um = self.ctx.z_scale_um

                astro = self.deconv._deconvolve_3d_stack(
                    astro, channel_label_key="channel.cells"
                )

                if self.ctx.has_dapi and dapi is not None and np.any(dapi > 0):
                    dapi = self.deconv._deconvolve_3d_stack_dapi(
                        dapi, pixel_xy_um, pixel_z_um,
                        channel_label_key="channel.dapi"
                    )

                if getattr(self.ctx, "has_granules", False) and granules_stack is not None:
                    granules_stack = self.deconv._deconvolve_3d_stack(
                        granules_stack, channel_label_key="channel.granules"
                    )

            gran_norm_clean = None
            if getattr(self.ctx, "has_granules", False) and granules_stack is not None:
                gran_norm_clean = self.segmenter._compute_astro_norm(
                    granules_stack, self.ctx.image_mode, is_granules=True
                )
                if gran_norm_clean is not None:
                    granules_stack = gran_norm_clean

            base = os.path.splitext(os.path.basename(filepath))[0]
            suffix = self.measurer._get_cell_suffix()
            out_dir = os.path.join(os.path.dirname(filepath), f"{base}_{suffix}_analysis")
            os.makedirs(out_dir, exist_ok=True)

            self.ctx._current_base_name = base
            self.ctx._current_out_dir = out_dir

            self.viz._save_deconvolved_stacks(astro, dapi, base, out_dir, granules=granules_stack)

            t2 = time.time()

            astro_norm_clean = self.segmenter._compute_astro_norm(astro, self.ctx.image_mode)
            astro_norm = astro_norm_clean

            if astro_norm is None:
                self.ctx.log(self.ctx.t("log.cells_no_signal"))

            if getattr(self.ctx, "auto_soma", False):
                soma_result = self.segmenter.auto_soma_setting(
                    astro, astro_norm=astro_norm
                )
                if soma_result is None:
                    self.ctx.log(self.ctx.t("log.soma_auto_none"))
                    soma_result = self.segmenter.manual_soma_setting(astro)
            else:
                soma_result = self.segmenter.manual_soma_setting(astro)

            self.ctx.soma_segmented = soma_result
            if self.ctx.soma_segmented.max() == 0:
                self.ctx.log(self.ctx.t("log.soma_not_marked"))
                return None, None

            click_centers = getattr(self.ctx, "_last_click_centers", []) or []
            click_types = getattr(self.ctx, "_last_click_types", []) or []

            if any(t == "donut" for t in click_types):
                self.ctx.soma_segmented = self.segmenter._detect_donut_and_fill(
                    astro_norm_clean, click_centers, click_types, self.ctx.soma_segmented
                )

                astro_norm_filled = astro_norm_clean.copy()
                astro_norm_filled[self.ctx.soma_segmented > 0] = 1.0
                astro_norm = astro_norm_filled

            if self.ctx._gui_do_filter:
                p0, p99 = np.percentile(astro_norm[astro_norm > 0], [1, 99])
                astro_norm_disp = np.clip((astro_norm - p0) / (p99 - p0 + 1e-6), 0, 1)
                mip_disp = np.max(astro_norm_disp, axis=0)

                soma_mip = np.max(self.ctx.soma_segmented, axis=0)
                unique_somas = np.unique(soma_mip)[1:]
                soma_colors = {}
                for i, cid in enumerate(unique_somas):
                    hue = (i * 0.618033988749895) % 1.0
                    r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                    soma_colors[int(cid)] = (r, g, b)

                while True:
                    res = self.ctx.bridge._ask(
                        "quality_check",
                        mip=mip_disp,
                        cell_mip=soma_mip,
                        colors=soma_colors,
                        mode="soma"
                    )

                    if res is None or res == "ok":
                        break

                    if res == "fused":
                        self.ctx.log(self.ctx.t("log.filter_fused_correction"))
                        current_centers = []
                        current_types = []
                        if self.ctx.soma_segmented is not None and self.ctx.soma_segmented.max() > 0:
                            soma_int = self.ctx.soma_segmented.astype(np.int32)
                            props = measure.regionprops(soma_int)
                            last_centers = getattr(self.ctx, "_last_click_centers", []) or []
                            last_types = getattr(self.ctx, "_last_click_types", []) or []
                            z_scale = float(self.ctx.z_scale_um) if self.ctx.scale_initialized else 1.0
                            xy_scale = float(self.ctx.xy_scale_um) * float(self.ctx.downsample) if self.ctx.scale_initialized else 1.0
                            for p in props:
                                pz, py, px = p.centroid
                                best = None
                                best_d = None
                                for lc, lt in zip(last_centers, last_types):
                                    lcz, lcy, lcx = lc
                                    dz = (lcz - pz) * z_scale
                                    dy = (lcy - py) * xy_scale
                                    dx = (lcx - px) * xy_scale
                                    d = dz * dz + dy * dy + dx * dx
                                    if best_d is None or d < best_d:
                                        best_d = d
                                        best = (int(lcy), int(lcx), int(lcz), lt)
                                if best is not None:
                                    current_centers.append((best[0], best[1], best[2]))
                                    current_types.append(best[3])
                                else:
                                    current_centers.append((int(py), int(px), int(pz)))
                                    current_types.append("solid")

                        astro_norm = astro_norm_clean.copy()
                        self.ctx.soma_segmented = np.zeros_like(self.ctx.soma_segmented, dtype=np.int32)

                        if getattr(self.ctx, "auto_soma", False):
                            new_soma = self.segmenter.auto_soma_setting(astro, astro_norm=astro_norm_clean, initial_points=current_centers, initial_types=current_types)
                            if new_soma is None:
                                new_soma = self.segmenter.manual_soma_setting(astro, initial_points=current_centers, initial_types=current_types)
                        else:
                            new_soma = self.segmenter.manual_soma_setting(astro, initial_points=current_centers, initial_types=current_types)

                        if new_soma is not None and new_soma.max() > 0:
                            self.ctx.soma_segmented = new_soma

                        click_centers = getattr(self.ctx, "_last_click_centers", []) or []
                        click_types = getattr(self.ctx, "_last_click_types", []) or []

                        if any(t == "donut" for t in click_types):
                            self.ctx.soma_segmented = self.segmenter._detect_donut_and_fill(
                                astro_norm_clean, click_centers, click_types, self.ctx.soma_segmented
                            )
                            astro_norm_filled = astro_norm_clean.copy()
                            astro_norm_filled[self.ctx.soma_segmented > 0] = 1.0
                            astro_norm = astro_norm_filled

                        p0, p99 = np.percentile(astro_norm[astro_norm > 0], [1, 99])
                        astro_norm_disp = np.clip((astro_norm - p0) / (p99 - p0 + 1e-6), 0, 1)
                        mip_disp = np.max(astro_norm_disp, axis=0)

                        soma_mip = np.max(self.ctx.soma_segmented, axis=0)
                        unique_somas = np.unique(soma_mip)[1:]
                        soma_colors = {}
                        for i, cid in enumerate(unique_somas):
                            hue = (i * 0.618033988749895) % 1.0
                            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                            soma_colors[int(cid)] = (r, g, b)
                        continue

            self.viz._save_slices_visualization(astro, self.ctx.soma_segmented,
                                                self.ctx.t("log.title_soma_slices"),
                                                "soma_slices")

            self.ctx.nuclei_segmented = self.segmenter.segment_nuclei_local_from_soma(
                self.ctx.soma_segmented, dapi
            )

            if self.ctx._gui_do_filter and self.ctx.has_dapi and self.ctx.nuclei_segmented.max() > 0:
                p0, p99 = np.percentile(dapi[dapi > 0], [1, 99])
                dapi_norm_disp = np.clip((dapi - p0) / (p99 - p0 + 1e-6), 0, 1)
                dapi_mip_disp = np.max(dapi_norm_disp, axis=0)

                nuclei_mip = np.max(self.ctx.nuclei_segmented, axis=0)
                unique_nuclei = np.unique(nuclei_mip)[1:]
                nuclei_colors = {}
                for i, cid in enumerate(unique_nuclei):
                    hue = (i * 0.618033988749895) % 1.0
                    r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                    nuclei_colors[int(cid)] = (r, g, b)

                while True:
                    res = self.ctx.bridge._ask(
                        "quality_check",
                        mip=dapi_mip_disp,
                        cell_mip=nuclei_mip,
                        colors=nuclei_colors,
                        mode="nuclei"
                    )

                    if res is None or res == "ok":
                        break

                    if res == "delete":
                        self.ctx.log(self.ctx.t("log.filter_manual_deletion"))
                        del_res = self.ctx.bridge._ask(
                            "cell_deletion",
                            mip=dapi_mip_disp,
                            cell_mip=nuclei_mip,
                            colors=nuclei_colors
                        )
                        if del_res:
                            for cid in del_res:
                                self.ctx.nuclei_segmented[self.ctx.nuclei_segmented == cid] = 0

                        nuclei_mip = np.max(self.ctx.nuclei_segmented, axis=0)
                        unique_nuclei = np.unique(nuclei_mip)[1:]
                        nuclei_colors = {}
                        for i, cid in enumerate(unique_nuclei):
                            hue = (i * 0.618033988749895) % 1.0
                            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                            nuclei_colors[int(cid)] = (r, g, b)
                        continue

                    if res == "artifacts":
                        if self.ctx.deconvolve_images:
                            self.ctx.log(self.ctx.t("log.filter_artifacts_correction"))
                            astro = self.deconv._deconvolve_3d_stack(
                                astro, channel_label_key="channel.cells"
                            )
                            if self.ctx.has_dapi and dapi is not None and np.any(dapi > 0):
                                dapi = self.deconv._deconvolve_3d_stack_dapi(
                                    dapi, pixel_xy_um, pixel_z_um,
                                    channel_label_key="channel.dapi"
                                )
                            if getattr(self.ctx, "has_granules", False) and granules_stack is not None:
                                granules_stack = self.deconv._deconvolve_3d_stack(
                                    granules_stack, channel_label_key="channel.granules"
                                )
                                gran_norm_clean = self.segmenter._compute_astro_norm(
                                    granules_stack, self.ctx.image_mode
                                )
                                if gran_norm_clean is not None:
                                    granules_stack = gran_norm_clean

                            astro_norm_clean = self.segmenter._compute_astro_norm(astro, self.ctx.image_mode)
                            astro_norm = astro_norm_clean

                            self.ctx.nuclei_segmented = self.segmenter.segment_nuclei_local_from_soma(self.ctx.soma_segmented, dapi)

                            p0, p99 = np.percentile(dapi[dapi > 0], [1, 99])
                            dapi_norm_disp = np.clip((dapi - p0) / (p99 - p0 + 1e-6), 0, 1)
                            dapi_mip_disp = np.max(dapi_norm_disp, axis=0)

                            nuclei_mip = np.max(self.ctx.nuclei_segmented, axis=0)
                            unique_nuclei = np.unique(nuclei_mip)[1:]
                            nuclei_colors = {}
                            for i, cid in enumerate(unique_nuclei):
                                hue = (i * 0.618033988749895) % 1.0
                                r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                                nuclei_colors[int(cid)] = (r, g, b)
                        else:
                            self.ctx.log(self.ctx.t("log.filter_artifacts_no_deconv"))
                            break
                        continue

            self.viz._save_slices_visualization(dapi, self.ctx.nuclei_segmented,
                                                self.ctx.t("log.title_nuclei_slices"),
                                                "nuclei_slices")

            if astro_norm_clean is not None:
                self.ctx.cells_segmented = self.segmenter._grow_cells_from_soma(
                    astro_norm_clean, self.ctx.soma_segmented
                )
            else:
                self.ctx.cells_segmented = np.zeros_like(self.ctx.soma_segmented,
                                                         dtype=np.int32)

            self.viz._save_slices_visualization(astro, self.ctx.cells_segmented,
                                                self.ctx.t("log.title_cells_slices"),
                                                "cells_slices")

            self.ctx.log(self.ctx.t("log.segmentation_time",
                                    t=time.time() - t2))

            self.ctx.bridge.app._deconv_was_on = self.ctx.deconvolve_images

            if self.ctx._gui_do_filter:
                p0, p99 = np.percentile(astro_norm[astro_norm > 0], [1, 99])
                astro_norm_disp = np.clip((astro_norm - p0) / (p99 - p0 + 1e-6), 0, 1)
                mip_disp = np.max(astro_norm_disp, axis=0)

                cell_mip = np.zeros((astro.shape[1], astro.shape[2]), dtype=np.int32)
                for z in range(astro.shape[0]):
                    cell_mip[self.ctx.cells_segmented[z] > 0] = self.ctx.cells_segmented[z][self.ctx.cells_segmented[z] > 0]

                unique_cells = np.unique(cell_mip)[1:]
                cell_colors = {}
                for i, cid in enumerate(unique_cells):
                    hue = (i * 0.618033988749895) % 1.0
                    r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                    cell_colors[int(cid)] = (r, g, b)

                while True:
                    res = self.ctx.bridge._ask(
                        "quality_check",
                        mip=mip_disp,
                        cell_mip=cell_mip,
                        colors=cell_colors,
                        mode="cells"
                    )

                    if res is None or res == "ok":
                        break

                    if res == "fused":
                        self.ctx.log(self.ctx.t("log.filter_fused_correction"))
                        current_centers = []
                        current_types = []
                        if self.ctx.soma_segmented is not None and self.ctx.soma_segmented.max() > 0:
                            props = measure.regionprops(self.ctx.soma_segmented)
                            last_centers = getattr(self.ctx, "_last_click_centers", []) or []
                            last_types = getattr(self.ctx, "_last_click_types", []) or []
                            z_scale = float(self.ctx.z_scale_um) if self.ctx.scale_initialized else 1.0
                            xy_scale = float(self.ctx.xy_scale_um) * float(self.ctx.downsample) if self.ctx.scale_initialized else 1.0
                            for p in props:
                                pz, py, px = p.centroid
                                best = None
                                best_d = None
                                for lc, lt in zip(last_centers, last_types):
                                    lcz, lcy, lcx = lc
                                    dz = (lcz - pz) * z_scale
                                    dy = (lcy - py) * xy_scale
                                    dx = (lcx - px) * xy_scale
                                    d = dz * dz + dy * dy + dx * dx
                                    if best_d is None or d < best_d:
                                        best_d = d
                                        best = (int(lcy), int(lcx), int(lcz), lt)
                                if best is not None:
                                    current_centers.append((best[0], best[1], best[2]))
                                    current_types.append(best[3])
                                else:
                                    current_centers.append((int(py), int(px), int(pz)))
                                    current_types.append("solid")

                        astro_norm = astro_norm_clean.copy()
                        self.ctx.soma_segmented = np.zeros_like(self.ctx.soma_segmented, dtype=np.int32)

                        if getattr(self.ctx, "auto_soma", False):
                            new_soma = self.segmenter.auto_soma_setting(astro, astro_norm=astro_norm_clean, initial_points=current_centers, initial_types=current_types)
                            if new_soma is None:
                                new_soma = self.segmenter.manual_soma_setting(astro, initial_points=current_centers, initial_types=current_types)
                        else:
                            new_soma = self.segmenter.manual_soma_setting(astro, initial_points=current_centers, initial_types=current_types)

                        if new_soma is not None and new_soma.max() > 0:
                            self.ctx.soma_segmented = new_soma
                            self.ctx.nuclei_segmented = self.segmenter.segment_nuclei_local_from_soma(self.ctx.soma_segmented, dapi)
                            self.ctx.cells_segmented = self.segmenter._grow_cells_from_soma(astro_norm_clean, self.ctx.soma_segmented)

                        click_centers = getattr(self.ctx, "_last_click_centers", []) or []
                        click_types = getattr(self.ctx, "_last_click_types", []) or []

                        if any(t == "donut" for t in click_types):
                            self.ctx.soma_segmented = self.segmenter._detect_donut_and_fill(
                                astro_norm_clean, click_centers, click_types, self.ctx.soma_segmented
                            )
                            astro_norm_filled = astro_norm_clean.copy()
                            astro_norm_filled[self.ctx.soma_segmented > 0] = 1.0
                            astro_norm = astro_norm_filled
                            self.ctx.cells_segmented = self.segmenter._grow_cells_from_soma(astro_norm, self.ctx.soma_segmented)

                        cell_mip = np.zeros((astro.shape[1], astro.shape[2]), dtype=np.int32)
                        for z in range(astro.shape[0]):
                            cell_mip[self.ctx.cells_segmented[z] > 0] = self.ctx.cells_segmented[z][self.ctx.cells_segmented[z] > 0]
                        unique_cells = np.unique(cell_mip)[1:]
                        cell_colors = {}
                        for i, cid in enumerate(unique_cells):
                            hue = (i * 0.618033988749895) % 1.0
                            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                            cell_colors[int(cid)] = (r, g, b)
                        continue

                    if res == "artifacts":
                        if self.ctx.deconvolve_images:
                            self.ctx.log(self.ctx.t("log.filter_artifacts_correction"))
                            astro = self.deconv._deconvolve_3d_stack(
                                astro, channel_label_key="channel.cells"
                            )
                            if self.ctx.has_dapi and dapi is not None and np.any(dapi > 0):
                                dapi = self.deconv._deconvolve_3d_stack_dapi(
                                    dapi, pixel_xy_um, pixel_z_um,
                                    channel_label_key="channel.dapi"
                                )
                            if getattr(self.ctx, "has_granules", False) and granules_stack is not None:
                                granules_stack = self.deconv._deconvolve_3d_stack(
                                    granules_stack, channel_label_key="channel.granules"
                                )
                                gran_norm_clean = self.segmenter._compute_astro_norm(
                                    granules_stack, self.ctx.image_mode
                                )
                                if gran_norm_clean is not None:
                                    granules_stack = gran_norm_clean

                            astro_norm_clean = self.segmenter._compute_astro_norm(astro, self.ctx.image_mode)
                            astro_norm = astro_norm_clean
                            if getattr(self.ctx, "auto_soma", False):
                                new_soma = self.segmenter.auto_soma_setting(astro, astro_norm=astro_norm_clean)
                                if new_soma is None:
                                    new_soma = self.segmenter.manual_soma_setting(astro)
                            else:
                                new_soma = self.segmenter.manual_soma_setting(astro)

                            if new_soma is not None and new_soma.max() > 0:
                                self.ctx.soma_segmented = new_soma
                                self.ctx.nuclei_segmented = self.segmenter.segment_nuclei_local_from_soma(self.ctx.soma_segmented, dapi)
                                self.ctx.cells_segmented = self.segmenter._grow_cells_from_soma(astro_norm_clean, self.ctx.soma_segmented)

                            p0, p99 = np.percentile(astro_norm[astro_norm > 0], [1, 99])
                            astro_norm_disp = np.clip((astro_norm - p0) / (p99 - p0 + 1e-6), 0, 1)
                            mip_disp = np.max(astro_norm_disp, axis=0)

                            cell_mip = np.zeros((astro.shape[1], astro.shape[2]), dtype=np.int32)
                            for z in range(astro.shape[0]):
                                cell_mip[self.ctx.cells_segmented[z] > 0] = self.ctx.cells_segmented[z][self.ctx.cells_segmented[z] > 0]
                            unique_cells = np.unique(cell_mip)[1:]
                            cell_colors = {}
                            for i, cid in enumerate(unique_cells):
                                hue = (i * 0.618033988749895) % 1.0
                                r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                                cell_colors[int(cid)] = (r, g, b)
                        else:
                            self.ctx.log(self.ctx.t("log.filter_artifacts_no_deconv"))
                            break
                        continue

                    if res == "cut":
                        self.ctx.log(self.ctx.t("log.filter_cut_soft_mode"))
                        self.ctx._growth_mode = "soft"
                        self.ctx.cells_segmented = self.segmenter._grow_cells_from_soma(astro_norm, self.ctx.soma_segmented)

                        cell_mip = np.zeros((astro.shape[1], astro.shape[2]), dtype=np.int32)
                        for z in range(astro.shape[0]):
                            cell_mip[self.ctx.cells_segmented[z] > 0] = self.ctx.cells_segmented[z][self.ctx.cells_segmented[z] > 0]
                        unique_cells = np.unique(cell_mip)[1:]
                        cell_colors = {}
                        for i, cid in enumerate(unique_cells):
                            hue = (i * 0.618033988749895) % 1.0
                            r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                            cell_colors[int(cid)] = (r, g, b)
                        continue

                    if res == "delete":
                        self.ctx.log(self.ctx.t("log.filter_manual_deletion"))
                        del_res = self.ctx.bridge._ask(
                            "cell_deletion",
                            mip=mip_disp,
                            cell_mip=cell_mip,
                            colors=cell_colors
                        )
                        if del_res:
                            for cid in del_res:
                                self.ctx.cells_segmented[self.ctx.cells_segmented == cid] = 0

                            cell_mip = np.zeros((astro.shape[1], astro.shape[2]), dtype=np.int32)
                            for z in range(astro.shape[0]):
                                cell_mip[self.ctx.cells_segmented[z] > 0] = self.ctx.cells_segmented[z][self.ctx.cells_segmented[z] > 0]

                            unique_cells = np.unique(cell_mip)[1:]
                            cell_colors = {}
                            for i, cid in enumerate(unique_cells):
                                hue = (i * 0.618033988749895) % 1.0
                                r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 1.0)
                                cell_colors[int(cid)] = (r, g, b)
                        continue

            self.ctx._growth_mode = "normal"

            gran_mask_all = None

            if self.ctx._gui_do_filter and granules_stack is not None:
                from scipy.ndimage import binary_opening

                gran_mask_all = np.zeros_like(self.ctx.cells_segmented, dtype=np.int32)
                unique_cells = np.unique(self.ctx.cells_segmented)[1:]

                for cid in unique_cells:
                    cell_mask = (self.ctx.cells_segmented == cid)
                    if not np.any(cell_mask):
                        continue

                    dilated_cell = binary_dilation(cell_mask, iterations=2)

                    gvals_dilated = granules_stack[dilated_cell]
                    if gvals_dilated.size < 10:
                        continue

                    gmax = float(np.max(gvals_dilated))
                    gmin = float(np.min(gvals_dilated))

                    if gmax > gmin:
                        try:
                            from skimage.filters import threshold_otsu
                            thr = float(threshold_otsu(gvals_dilated))
                        except Exception:
                            thr = gmin + (gmax - gmin) * 0.5

                        thr = max(thr, float(np.percentile(gvals_dilated, 50)))

                        c_gran_mask = dilated_cell & (granules_stack >= thr)
                        c_gran_mask = binary_opening(c_gran_mask, structure=np.ones((2, 2, 2)))

                        lbl_g, n_g = label(c_gran_mask, structure=np.ones((3, 3, 3)))

                        if n_g > 0:
                            sizes = np.bincount(lbl_g.ravel())
                            valid_labels = np.where((sizes >= 5) & (sizes <= 2000))[0]
                            valid_labels = valid_labels[valid_labels > 0]

                            for lid in valid_labels:
                                granule_coords = np.argwhere(lbl_g == lid)
                                if len(granule_coords) == 0:
                                    continue

                                center_idx = len(granule_coords) // 2
                                cz, cy, cx = granule_coords[center_idx]

                                if cell_mask[cz, cy, cx]:
                                    gran_mask_all[lbl_g == lid] = cid

                gran_mip = np.max(gran_mask_all, axis=0)
                gran_colors = cell_colors

                while True:
                    res = self.ctx.bridge._ask(
                        "quality_check",
                        mip=mip_disp,
                        cell_mip=cell_mip,
                        gran_mip=gran_mip,
                        colors=gran_colors,
                        mode="granules"
                    )

                    if res is None or res == "ok":
                        break

                    if res == "delete":
                        self.ctx.log(self.ctx.t("log.filter_manual_deletion"))
                        del_res = self.ctx.bridge._ask(
                            "cell_deletion",
                            mip=mip_disp,
                            cell_mip=gran_mip,
                            colors=gran_colors
                        )
                        if del_res:
                            for cid in del_res:
                                gran_mask_all[gran_mask_all == cid] = 0

                        gran_mip = np.max(gran_mask_all, axis=0)
                        continue

            self.viz._save_cells_overlay_stack(astro, self.ctx.cells_segmented, base, out_dir)

            self.viz.save_3d_stacks(
                astro, dapi, base, out_dir,
                soma=self.ctx.soma_segmented,
                nuclei=self.ctx.nuclei_segmented,
                cells=self.ctx.cells_segmented,
            )

            t3 = time.time()
            measurements = self.measurer.measure_cells_3d(
                self.ctx.cells_segmented, astro_norm_clean, self.ctx.nuclei_segmented,
                filepath,
                granules=gran_mask_all,
                cell_channel=astro,
            )
            self.ctx.log(self.ctx.t("log.measurement_time",
                                    t=time.time() - t3))

            self.viz.save_visualizations(astro, dapi, self.ctx.cells_segmented,
                                         self.ctx.nuclei_segmented, base, out_dir)

            if len(measurements) > 0:
                measurements.to_csv(
                    os.path.join(out_dir, f"{base}_measurements.csv"),
                    index=False
                )
                self.measurer._append_to_all_measurements(measurements, filepath)

            self.measurer._append_summary_count(filepath, len(measurements))

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
            if kind == "channel_selector": self._dialog_channel_selector(kwargs)
            elif kind == "soma_picker":      self._dialog_soma_picker(kwargs)
            elif kind == "soma_confirm":     self._dialog_soma_confirm(kwargs)
            elif kind == "slice_selector":   self._dialog_slice_selector(kwargs)
            elif kind == "quality_check":    self._dialog_quality_check(kwargs)
            elif kind == "cell_deletion":    self._dialog_cell_deletion(kwargs)
            elif kind == "preprocess_choice": self._dialog_preprocess_choice(kwargs)
            elif kind == "crop_editor":       self._dialog_crop_editor(kwargs)
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

    def _dialog_metrics_selector(self, cfg_ref):
        panel = self.app.start_work_panel(self.app._t("dlg.metrics.title"))

        tk.Label(panel, text=self.app._t("dlg.metrics.hint"),
                 bg=self.app.BG_CARD, fg=self.app.FG_DIM,
                 font=("Segoe UI", 9, "italic"),
                 justify="left", anchor="w", wraplength=480).pack(
            anchor="w", pady=(0, 8))

        top_bar = tk.Frame(panel, bg=self.app.BG_CARD)
        top_bar.pack(fill="x", pady=(0, 6))

        vars_map = {}

        for _grp, items in METRIC_GROUPS:
            for it in items:
                vars_map[it] = tk.BooleanVar(value=bool(cfg_ref.get(it, True)))

        def _set_all(val):
            for v in vars_map.values():
                v.set(bool(val))
            self.app._pending_metrics_config = {k: bool(v.get()) for k, v in vars_map.items()}

        RoundedButton(top_bar, self.app._t("dlg.metrics.select_all"),
                      command=lambda: _set_all(True),
                      width=130, height=28, radius=14,
                      bg=self.app.BG_HOVER, fg=self.app.FG,
                      hover_bg=self.app.ACCENT, hover_fg="#1a1c22",
                      active_bg=self.app.ACCENT_D,
                      font=("Segoe UI", 9)).pack(side="left", padx=(0, 6))
        RoundedButton(top_bar, self.app._t("dlg.metrics.select_none"),
                      command=lambda: _set_all(False),
                      width=130, height=28, radius=14,
                      bg=self.app.BG_HOVER, fg=self.app.FG,
                      hover_bg=self.app.ACCENT, hover_fg="#1a1c22",
                      active_bg=self.app.ACCENT_D,
                      font=("Segoe UI", 9)).pack(side="left")

        scroll_wrap = tk.Frame(panel, bg=self.app.BG_CARD)
        scroll_wrap.pack(fill="both", expand=True)

        canvas = tk.Canvas(scroll_wrap, bg=self.app.BG_CARD, highlightthickness=0)
        sb = RoundedScrollbar(scroll_wrap, command=canvas.yview,
                              width=10, orient="vertical",
                              bg=self.app.BG_CARD,
                              thumb=self.app.BORDER,
                              thumb_hover=self.app.ACCENT)
        canvas.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        inner = tk.Frame(canvas, bg=self.app.BG_CARD)
        canvas.create_window((0, 0), window=inner, anchor="nw")

        def _on_frame_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))

        inner.bind("<Configure>", _on_frame_configure)

        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        def _on_var_change(*_):
            self.app._pending_metrics_config = {k: bool(v.get()) for k, v in vars_map.items()}

        is_2d = bool(getattr(self.app, "_probe_is_2d", False))

        for grp_key, items in METRIC_GROUPS:
            if grp_key.endswith("_2d_cell") or grp_key.endswith("_2d_soma") or grp_key.endswith("_2d_nucleus"):
                grp_is_2d = True
            elif grp_key.endswith("_3d_cell") or grp_key.endswith("_3d_soma") or grp_key.endswith("_3d_nucleus"):
                grp_is_2d = False
            else:
                grp_is_2d = None

            if grp_is_2d is not None:
                if is_2d and not grp_is_2d:
                    continue
                if (not is_2d) and grp_is_2d:
                    continue

            grp_frame = tk.LabelFrame(
                inner, text=self.app._t("dlg.metrics." + grp_key),
                bg=self.app.BG_CARD, fg=self.app.ACCENT,
                font=("Segoe UI", 9, "bold"),
                bd=0, highlightthickness=1,
                highlightbackground=self.app.BORDER,
            )
            grp_frame.pack(fill="x", pady=4, padx=2)

            for it in items:
                key = "dlg.metrics.item." + it
                chk = tk.Checkbutton(
                    grp_frame,
                    text=self.app._t(key),
                    variable=vars_map[it],
                    command=_on_var_change,
                    bg=self.app.BG_CARD, fg=self.app.FG,
                    activebackground=self.app.BG_CARD,
                    activeforeground=self.app.ACCENT,
                    selectcolor=self.app.BG_SOFT,
                    relief="flat", bd=0,
                    font=("Segoe UI", 9),
                    anchor="w", justify="left",
                    highlightthickness=0,
                    cursor="hand2",
                )
                chk.pack(anchor="w", padx=8, pady=1)

        def finish():
            canvas.unbind_all("<MouseWheel>")
            for k, v in vars_map.items():
                cfg_ref[k] = bool(v.get())
            self.app._pending_metrics_config = dict(cfg_ref)

            ch_cfg = getattr(self.app, "_pending_channels_config", None)
            if ch_cfg:
                ch_cfg["metrics_config"] = dict(cfg_ref)

            self.app.stop_work_panel()

            if getattr(self.app, "_active_channel_finish", None) is not None:
                self.app._active_channel_finish()
                self.app._active_channel_finish = None
            else:
                self._ok(None)

        def cancel():
            canvas.unbind_all("<MouseWheel>")
            self.app.stop_work_panel()
            self._ok(None)

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            cancel=True, cancel_value=None,
            cancel_text=self.app._t("btn.cancel"),
            help_text="",
            on_cancel=cancel,
        )

    def _dialog_channel_selector(self, kw):
        mips = kw.get("mips", [])
        n = len(mips)
        if n < 1:
            self._ok({
                "selected_channels": [0],
                "cell_idx": 0,
                "has_dapi": False,
                "dapi_idx": None,
                "has_granules": False,
                "granule_idx": None,
                "metrics_config": default_metrics_config(),
            })
            return

        filepath = getattr(self.app, "_current_filepath", None)
        cell_default = kw.get("cell_idx_default", 0) or 0
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

        scroll_frame = tk.Frame(panel, bg=self.app.BG_CARD)
        scroll_frame.pack(fill="both", expand=True)

        canvas = tk.Canvas(scroll_frame, bg=self.app.BG_CARD, highlightthickness=0)
        scrollbar = RoundedScrollbar(scroll_frame, command=canvas.xview,
                                     width=10, orient="horizontal",
                                     bg=self.app.BG_CARD,
                                     thumb=self.app.BORDER,
                                     thumb_hover=self.app.ACCENT)
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

        roles = ["none"] * n
        roles[int(cell_default)] = "cells"
        if has_dapi_default and dapi_default is not None:
            roles[int(dapi_default)] = "dapi"

        role_vars = [tk.StringVar(value=roles[i]) for i in range(n)]
        tile_drawers = {}
        dropdowns = {}

        role_to_key = {
            "cells": "dlg.channels.role_cells",
            "dapi": "dlg.channels.role_dapi",
            "granules": "dlg.channels.role_granules",
            "none": "dlg.channels.role_none",
        }
        role_to_color = {
            "cells": self.app.OK,
            "dapi": self.app.ACCENT,
            "granules": self.app.WARN,
            "none": self.app.FG_MUTED,
        }

        def _roles_in_use():
            used = {"cells": None, "dapi": None, "granules": None}
            for i in range(n):
                r = role_vars[i].get()
                if r in used:
                    used[r] = i
            return used

        def _refresh_all():
            used = _roles_in_use()
            for i in range(n):
                if i in tile_drawers:
                    tile_drawers[i]()
                if i in dropdowns:
                    dropdowns[i].set_state(role_vars[i].get(), used, i)

        def _update_pending_config():
            used = _roles_in_use()
            cell_idx = used.get("cells")
            if cell_idx is not None:
                dapi_idx = used.get("dapi")
                granule_idx = used.get("granules")
                selected = sorted({i for i in range(n) if role_vars[i].get() != "none"} | {cell_idx})
                self.app._pending_channels_config = {
                    "selected_channels": selected,
                    "cell_idx": int(cell_idx),
                    "has_dapi": bool(dapi_idx is not None),
                    "dapi_idx": int(dapi_idx) if dapi_idx is not None else None,
                    "has_granules": bool(granule_idx is not None),
                    "granule_idx": int(granule_idx) if granule_idx is not None else None,
                    "metrics_config": dict(self.app._pending_metrics_config),
                }

        def _set_role(index, new_role):
            old = role_vars[index].get()
            if new_role == old:
                return

            if new_role != "none":
                for k in range(n):
                    if k != index and role_vars[k].get() == new_role:
                        role_vars[k].set("none")
                        break

            if new_role == "none" and old == "cells":
                for k in range(n):
                    if k != index:
                        role_vars[k].set("cells")
                        break

            role_vars[index].set(new_role)
            _refresh_all()
            _update_pending_config()

            used = _roles_in_use()
            if filepath:
                self.app._update_optics_for_channels(filepath, used.get("cells"), used.get("dapi"), used.get("granules"))

        def _make_tile_column(index, mip):
            arr = self._to_uint8(mip)
            pil = self._pil_from_gray(arr)
            w, h, _ = self._fit_size(pil.width, pil.height, max_w=200, max_h=200)
            pil = pil.resize((w, h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)

            column = tk.Frame(imgs_frame, bg=self.app.BG_CARD)
            column.pack(side="left", padx=6, pady=2, anchor="n")

            tile_w = w + 8
            tile_h = h + 8 + 56

            tile = tk.Canvas(column, width=tile_w, height=tile_h,
                             bg=self.app.BG_CARD,
                             highlightthickness=0, bd=0)
            tile.pack(side="top")
            tile.image = photo

            def draw_tile():
                tile.delete("all")
                r = role_vars[index].get()
                border = role_to_color.get(r, self.app.FG_MUTED)
                if r == "none":
                    border = self.app.BORDER

                tile.create_rectangle(1, 1, tile_w - 1, tile_h - 1,
                                      outline=border, width=2,
                                      fill=self.app.BG_SOFT)
                tile.create_image(4, 4, anchor="nw", image=photo)

                tile.create_text(tile_w // 2, tile_h - 44,
                                 text=self.app._t("dlg.channels.channel_n", n=index),
                                 fill=self.app.FG,
                                 font=("Segoe UI", 10, "bold"))
                role_lbl = self.app._t(role_to_key.get(r, "dlg.channels.role_none"))
                tile.create_text(tile_w // 2, tile_h - 22,
                                 text=role_lbl,
                                 fill=role_to_color.get(r, self.app.FG_DIM),
                                 font=("Segoe UI", 9, "bold"))

            tile_drawers[index] = draw_tile
            draw_tile()

            dropdown = RoleDropdown(
                column,
                app=self.app,
                on_change=lambda new_role, idx=index: _set_role(idx, new_role),
                labels={
                    "cells": self.app._t("dlg.channels.role_cells"),
                    "dapi": self.app._t("dlg.channels.role_dapi"),
                    "granules": self.app._t("dlg.channels.role_granules"),
                    "none": self.app._t("dlg.channels.role_none"),
                },
                taken_label=self.app._t("dlg.channels.role_used"),
                width=tile_w - 8,
            )
            dropdown.pack(side="top", pady=(4, 0))
            dropdowns[index] = dropdown

        for i in range(n):
            _make_tile_column(i, mips[i])

        _refresh_all()
        _update_pending_config()

        used = _roles_in_use()
        if filepath:
            self.app._update_optics_for_channels(filepath, used.get("cells"), used.get("dapi"), used.get("granules"))

        metrics_bar = tk.Frame(panel, bg=self.app.BG_CARD)
        metrics_bar.pack(anchor="w", pady=(8, 0))

        def _open_metrics_dialog():
            self.app._probe_is_2d = bool(getattr(self.app, "_probe_is_2d", False))
            self._dialog_metrics_selector(self.app._pending_metrics_config)

        RoundedButton(
            metrics_bar,
            self.app._t("dlg.metrics.btn"),
            command=_open_metrics_dialog,
            width=160, height=30, radius=15,
            bg=self.app.BG_HOVER, fg=self.app.FG,
            hover_bg=self.app.ACCENT, hover_fg="#1a1c22",
            active_bg=self.app.ACCENT_D,
            font=("Segoe UI", 9),
        ).pack(side="left")

        def finish():
            canvas.unbind_all("<MouseWheel>")
            cfg = getattr(self.app, "_pending_channels_config", None)
            if not cfg or cfg.get("cell_idx") is None:
                messagebox.showwarning(self.app._t("dlg.channels.title"),
                                       self.app._t("dlg.channels.choose_one"))
                return

            self.app._pending_metrics_config = dict(self.app._pending_metrics_config)
            if "metrics_config" in cfg:
                cfg["metrics_config"] = dict(self.app._pending_metrics_config)

            self.app.stop_work_panel()
            self._ok(cfg)

        def cancel():
            canvas.unbind_all("<MouseWheel>")
            self.app.stop_work_panel()
            self._ok(None)

        self.app._active_channel_finish = finish

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            cancel=True, cancel_value=None,
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.channels.help"),
            on_cancel=cancel,
        )

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
        click_types = kw.get("click_types", ["solid"] * len(click_centers))

        panel = self.app.start_work_panel(self.app._t("dlg.soma_confirm.title"))

        z_dim, img_h, img_w = astro_norm.shape
        z_var = tk.IntVar(value=0)
        photo_ref = {"photo": None}

        to_manual_bar = tk.Frame(panel, bg=self.app.BG_CARD)
        to_manual_bar.pack(side="bottom", fill="x", pady=(0, 6))

        def to_manual():
            self.app.stop_work_panel()
            self._ok("to_manual")

        RoundedButton(
            to_manual_bar,
            self.app._t("btn.to_manual_soma"),
            command=to_manual,
            width=320, height=32, radius=16,
            bg=self.app.BG_SOFT, fg=self.app.FG,
            hover_bg=self.app.ACCENT, hover_fg="#1a1c22",
            active_bg=self.app.ACCENT_D,
            font=("Segoe UI", 9)
        ).pack(side="left")

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
            for i, (cz, cy, cx) in enumerate(click_centers):
                if cz == z:
                    n_here += 1
                    px = x0 + cx * scale_x
                    py = y0 + cy * scale_y
                    r_ = 6
                    color = "#8fd694" if click_types[i] == "solid" else "#ffd54a"
                    canvas.create_oval(px - r_, py - r_, px + r_, py + r_,
                                       outline="white", fill=color, width=2)

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
        initial_points = kw.get("initial_points", None)
        initial_types = kw.get("initial_types", None)

        panel = self.app.start_work_panel(self.app._t("dlg.soma.title"))

        if astro_norm is not None:
            color_mip = self._build_color_mip(astro_norm)
        else:
            u8 = self._to_uint8(mip)
            color_mip = np.stack([u8, u8, u8], axis=-1)

        img_h, img_w = color_mip.shape[:2]

        points = []
        point_types = []
        if initial_points:
            z_max_allowed = astro_norm.shape[0] - 1 if astro_norm is not None else 0
            for i, pt in enumerate(initial_points):
                if isinstance(pt, (tuple, list)) and len(pt) == 3:
                    py_, px_, pz_ = pt
                    py_ = max(0, min(img_h - 1, int(py_)))
                    px_ = max(0, min(img_w - 1, int(px_)))
                    pz_ = max(0, min(z_max_allowed, int(pz_)))
                    points.append((py_, px_, pz_))
                    if initial_types and i < len(initial_types):
                        point_types.append(initial_types[i])
                    else:
                        point_types.append("solid")

        scale_ref = {"s": 1.0, "x0": 0, "y0": 0}
        photo_ref = {"photo": None}
        cursor_ref = {"x": None, "y": None}

        var_new_type = tk.StringVar(value="solid")

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

        counter = tk.Label(panel, text=self.app._t("dlg.soma.counter", n=len(points)),
                           bg=self.app.BG_CARD, fg=self.app.FG_DIM,
                           font=("Segoe UI", 9))
        counter.pack(side="bottom", anchor="w", pady=(6, 0))

        type_bar = tk.Frame(panel, bg=self.app.BG_CARD)
        type_bar.pack(side="bottom", fill="x", pady=(4, 0))
        tk.Label(type_bar, text=self.app._t("dlg.soma.type_label"), bg=self.app.BG_CARD, fg=self.app.FG, font=("Segoe UI", 9)).pack(side="left")
        tk.Radiobutton(type_bar, text=self.app._t("dlg.soma.type_solid"), variable=var_new_type, value="solid",
                       bg=self.app.BG_CARD, fg=self.app.FG, selectcolor=self.app.BG_SOFT, activebackground=self.app.BG_CARD,
                       font=("Segoe UI", 9)).pack(side="left", padx=4)
        tk.Radiobutton(type_bar, text=self.app._t("dlg.soma.type_donut"), variable=var_new_type, value="donut",
                       bg=self.app.BG_CARD, fg=self.app.FG, selectcolor=self.app.BG_SOFT, activebackground=self.app.BG_CARD,
                       font=("Segoe UI", 9)).pack(side="left", padx=4)

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

            canvas.create_line(0, y, cw, y, fill="#ffd54a", width=1, tags="crosshair")
            canvas.create_line(x, 0, x, ch, fill="#ffd54a", width=1, tags="crosshair")
            r = 3
            canvas.create_oval(x - r, y - r, x + r, y + r, outline="#1a1c22", fill="#ffd54a", width=1, tags="crosshair")

            coords_text = _get_cursor_orig_text(x, y)
            if coords_text is not None:
                canvas.create_text(10, 10, anchor="nw", text=coords_text, fill="#ffd54a", font=("Consolas", 10, "bold"), tags="crosshair")

        def _get_cursor_orig_text(sx, sy):
            s = scale_ref["s"]
            x0 = scale_ref["x0"]
            y0 = scale_ref["y0"]
            if s <= 0: return None
            ox = int((sx - x0) / s)
            oy = int((sy - y0) / s)
            if not (0 <= ox < img_w and 0 <= oy < img_h): return None
            return self.app._t("dlg.soma.coords", y=oy, x=ox)

        def _draw_overlay():
            if not overlay_ref["active"]: return
            ov_canvas = overlay_ref["canvas"]
            if ov_canvas is None: return

            crop = overlay_ref["crop"]
            z_var = overlay_ref["z_var"]
            z_lbl = overlay_ref["z_lbl"]
            crop_origin = overlay_ref["crop_origin"]
            crop_size = overlay_ref["crop_size"]
            if (crop is None or z_var is None or crop_origin is None or crop_size is None): return

            z_cur = int(z_var.get())
            z_dim_crop = crop.shape[0]
            z_cur = max(0, min(z_dim_crop - 1, z_cur))

            sl = crop[z_cur]
            nz = sl[sl > 0]
            if nz.size > 0:
                p2, p98 = np.percentile(nz, [2, 98])
                if p98 - p2 < 1e-6: sl_n = np.zeros_like(sl)
                else: sl_n = np.clip((sl - p2) / (p98 - p2), 0.0, 1.0)
            else:
                sl_n = np.zeros_like(sl)

            u8 = (sl_n * 255).astype(np.uint8)
            pil = self._pil_from_gray(u8).resize((OV_IMG, OV_IMG), Image.BILINEAR)
            base_rgb = np.array(pil.convert("RGB"))

            origin_y, origin_x = crop_origin
            half = crop_size // 2

            for i, pt in enumerate(points):
                py, px = pt[0], pt[1]
                dy = py - origin_y + half
                dx = px - origin_x + half
                if 0 <= dy < crop_size and 0 <= dx < crop_size:
                    cy_ = int(dy * OV_IMG / crop_size)
                    cx_ = int(dx * OV_IMG / crop_size)
                    r_ = 4
                    y0_ = max(0, cy_ - r_); y1_ = min(OV_IMG, cy_ + r_)
                    x0_ = max(0, cx_ - r_); x1_ = min(OV_IMG, cx_ + r_)
                    if y1_ > y0_ and x1_ > x0_:
                        color = [143, 214, 148] if point_types[i] == "solid" else [255, 213, 74]
                        base_rgb[y0_:y1_, x0_:x1_, 0] = color[0]
                        base_rgb[y0_:y1_, x0_:x1_, 1] = color[1]
                        base_rgb[y0_:y1_, x0_:x1_, 2] = color[2]

            photo = ImageTk.PhotoImage(Image.fromarray(base_rgb, mode="RGB"))
            overlay_ref["photo"] = photo

            ov_canvas.delete("all")
            ov_canvas.create_image(0, 0, anchor="nw", image=photo)

            if z_lbl is not None:
                z_lbl.config(text=self.app._t("dlg.soma.overlay_z", z=z_cur))

        def _close_overlay():
            if not overlay_ref["active"]: return
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
            if not overlay_ref["active"]: return
            z_var = overlay_ref["z_var"]
            crop_origin = overlay_ref["crop_origin"]
            crop_size = overlay_ref["crop_size"]
            if (z_var is None or crop_origin is None or crop_size is None): return
            if not (0 <= event.x < OV_IMG and 0 <= event.y < OV_IMG): return

            cx_in_crop = int(event.x * crop_size / OV_IMG)
            cy_in_crop = int(event.y * crop_size / OV_IMG)
            origin_y, origin_x = crop_origin
            half = crop_size // 2
            orig_y = origin_y - half + cy_in_crop
            orig_x = origin_x - half + cx_in_crop
            if not (0 <= orig_y < img_h and 0 <= orig_x < img_w): return

            z_cur = int(z_var.get())

            replaced = False
            for idx, pt in enumerate(points):
                if pt[0] == orig_y and pt[1] == orig_x:
                    points[idx] = (orig_y, orig_x, z_cur)
                    point_types[idx] = var_new_type.get()
                    replaced = True
                    break
            if not replaced:
                points.append((orig_y, orig_x, z_cur))
                point_types.append(var_new_type.get())

            counter.config(text=self.app._t("dlg.soma.counter", n=len(points)))
            redraw()
            _draw_overlay()

        def redraw():
            cw = canvas.winfo_width()
            ch = canvas.winfo_height()
            if cw <= 2 or ch <= 2: return

            s = min(cw / img_w, ch / img_h)
            new_w = max(1, int(img_w * s))
            new_h = max(1, int(img_h * s))
            x0 = (cw - new_w) // 2
            y0 = (ch - new_h) // 2

            scale_ref["s"] = s
            scale_ref["x0"] = x0
            scale_ref["y0"] = y0

            pil = Image.fromarray(color_mip, mode="RGB").resize((new_w, new_h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)
            photo_ref["photo"] = photo

            canvas.delete("all")
            canvas.create_image(x0, y0, anchor="nw", image=photo)

            r = 6
            for i, pt in enumerate(points):
                py, px = pt[0], pt[1]
                dx = x0 + px * s
                dy = y0 + py * s
                color = "#8fd694" if point_types[i] == "solid" else "#ffd54a"
                canvas.create_oval(dx - r, dy - r, dx + r, dy + r, outline="white", fill=color, width=2)

            _draw_crosshair()

        canvas.bind("<Configure>", lambda e: redraw())

        def _get_orig_coords(event):
            s = scale_ref["s"]
            x0 = scale_ref["x0"]
            y0 = scale_ref["y0"]
            if s <= 0: return None
            x = int((event.x - x0) / s)
            y = int((event.y - y0) / s)
            if not (0 <= x < img_w and 0 <= y < img_h): return None
            return y, x

        def _auto_pick_z(orig_y, orig_x):
            z_cur = 0
            best_score = -1.0
            if astro_norm is not None:
                for zz in range(astro_norm.shape[0]):
                    win = astro_norm[zz, max(0, orig_y - 2):min(img_h, orig_y + 3), max(0, orig_x - 2):min(img_w, orig_x + 3)]
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
            if overlay_ref["active"]: return
            coords = _get_orig_coords(event)
            if coords is None: return
            y, x = coords

            shift_held = bool(event.state & 0x0001)
            ctrl_held = bool(event.state & 0x0004)

            if ctrl_held:
                if points:
                    d = [(pt[0] - y) ** 2 + (pt[1] - x) ** 2 for pt in points]
                    idx = int(np.argmin(d))
                    points.pop(idx)
                    point_types.pop(idx)
            elif shift_held:
                if points:
                    d = [(pt[0] - y) ** 2 + (pt[1] - x) ** 2 for pt in points]
                    idx = int(np.argmin(d))
                    point_types[idx] = "donut" if point_types[idx] == "solid" else "solid"
            else:
                z_cur = _auto_pick_z(y, x)
                points.append((y, x, z_cur))
                point_types.append(var_new_type.get())

            counter.config(text=self.app._t("dlg.soma.counter", n=len(points)))
            redraw()

        def on_right(event):
            if overlay_ref["active"]: return
            coords = _get_orig_coords(event)
            if coords is None: return
            y, x = coords
            _open_overlay(y, x)

        canvas.bind("<Button-1>", on_left)
        canvas.bind("<Button-3>", on_right)
        canvas.bind("<Motion>", on_motion)
        canvas.bind("<Leave>", on_leave)

        def _open_overlay(orig_y, orig_x):
            _close_overlay()
            if astro_norm is None: return

            w = max(40, min(img_h, img_w) // 10)
            side = 2 * w
            crop = np.zeros((astro_norm.shape[0], side, side), dtype=astro_norm.dtype)

            y0c, y1c = orig_y - w, orig_y + w
            x0c, x1c = orig_x - w, orig_x + w

            sy0 = max(0, -y0c); sy1 = sy0 + (min(y1c, img_h) - max(y0c, 0))
            sx0 = max(0, -x0c); sx1 = sx0 + (min(x1c, img_w) - max(x0c, 0))
            dy0 = max(0, y0c); dy1 = dy0 + (sy1 - sy0)
            dx0 = max(0, x0c); dx1 = dx0 + (sx1 - sx0)

            crop[:, sy0:sy1, sx0:sx1] = astro_norm[:, dy0:dy1, dx0:dx1]

            crop_size = side
            crop_origin = (orig_y, orig_x)

            ov_frame = tk.Frame(canvas, bg=self.app.BG_CARD, highlightthickness=2, highlightbackground=self.app.ACCENT)
            ov_canvas = tk.Canvas(ov_frame, bg="#101216", highlightthickness=0, bd=0, width=OV_IMG, height=OV_IMG)
            ov_canvas.pack(side="top", padx=0, pady=0)

            z_max = crop.shape[0] - 1
            z_best = 0
            best_score = -1.0
            for zz in range(crop.shape[0]):
                sl = crop[zz]
                if sl.size == 0: continue
                score = float(np.max(sl))
                if score > best_score:
                    best_score = score
                    z_best = zz

            z_var = tk.IntVar(value=int(z_best))
            slider_holder = tk.Frame(ov_frame, bg=self.app.BG_CARD)
            slider_holder.pack(side="top", fill="x", padx=OV_PAD, pady=(4, 0))

            z_lbl = tk.Label(slider_holder, text=self.app._t("dlg.soma.overlay_z", z=int(z_best)), bg=self.app.BG_CARD, fg=self.app.FG, font=("Consolas", 9, "bold"))
            z_lbl.pack(side="top", anchor="w")

            RoundedSlider(slider_holder, from_=0, to=max(0, z_max), variable=z_var, length=OV_IMG - 20, command=_on_overlay_slide,
                          bg=self.app.BG_CARD, track=self.app.BORDER, fill=self.app.ACCENT, thumb="#d4d4dc", thumb_hover="#ffffff", fg=self.app.FG, showvalue=False).pack(fill="x", pady=(2, 0))

            btn_holder = tk.Frame(ov_frame, bg=self.app.BG_CARD)
            btn_holder.pack(side="bottom", fill="x", padx=OV_PAD, pady=(4, 6))
            RoundedButton(btn_holder, self.app._t("btn.ok"), command=_close_overlay, width=80, height=26, radius=13,
                          bg=self.app.ACCENT, fg="#1a1c22", hover_bg=self.app.ACCENT_D, hover_fg="#ffffff", active_bg="#4a6fbf", font=("Segoe UI", 9, "bold")).pack(side="right")

            ov_frame.update_idletasks()
            need_w = ov_frame.winfo_reqwidth()
            need_h = ov_frame.winfo_reqheight()
            cw = canvas.winfo_width()

            ov_x = cw - need_w - 10
            ov_y = 10
            if ov_x < 0: ov_x = 0

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
            self._ok((list(points), list(point_types)))

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

    def _dialog_quality_check(self, kw):
        mip = kw["mip"]
        cell_mip = kw["cell_mip"]
        colors = kw["colors"]
        mode = kw.get("mode", "cells")
        gran_mip = kw.get("gran_mip", None)

        if mode == "soma":
            entity = self.app._t("entity.soma")
        elif mode == "nuclei":
            entity = self.app._t("entity.nuclei")
        elif mode == "granules":
            entity = self.app._t("entity.granules")
        else:
            entity = self.app._t("entity.cells")

        y_dim, x_dim = mip.shape
        panel = self.app.start_work_panel(self.app._t("dlg.filter.quality_title"))

        question_text = self.app._t("dlg.filter.quality_ask_entity", entity=entity)
        tk.Label(panel, text=question_text,
                 bg=self.app.BG_CARD, fg=self.app.FG,
                 font=("Segoe UI", 10, "bold"), justify="left", anchor="w").pack(anchor="w", pady=(0, 10))

        base_u8 = self._to_uint8(mip)
        disp_w, disp_h, scale = self._fit_size(x_dim, y_dim, max_w=500, max_h=470)
        base_small = np.array(
            self._pil_from_gray(base_u8).resize((disp_w, disp_h), Image.BILINEAR))

        canvas = tk.Canvas(panel, width=disp_w, height=disp_h, bg="black",
                           highlightthickness=0)
        canvas.pack(pady=4)

        masks_small = {}
        for cid in np.unique(cell_mip)[1:]:
            m = (cell_mip == cid)
            if not np.any(m):
                continue
            small = np.array(
                Image.fromarray(m.astype(np.uint8) * 255, mode="L")
                .resize((disp_w, disp_h), Image.NEAREST)) > 127
            masks_small[cid] = small

        rgb = np.stack([base_small, base_small, base_small], axis=-1).astype(np.float32)
        for cid, small in masks_small.items():
            if not small.any():
                continue
            col = np.array(colors.get(cid, (1.0, 0.5, 0.0)), dtype=np.float32) * 255.0
            a = 0.25 if mode == "granules" else 0.55
            rgb[small] = rgb[small] * (1 - a) + col * a

        if mode == "granules" and gran_mip is not None:
            for cid in np.unique(gran_mip)[1:]:
                m_g = (gran_mip == cid)
                if not np.any(m_g):
                    continue
                small_g = np.array(
                    Image.fromarray(m_g.astype(np.uint8) * 255, mode="L")
                    .resize((disp_w, disp_h), Image.NEAREST)) > 127

                col = np.array(colors.get(cid, (1.0, 0.5, 0.0)), dtype=np.float32) * 255.0
                a = 0.95
                rgb[small_g] = rgb[small_g] * (1 - a) + col * a

        pil = Image.fromarray(rgb.astype(np.uint8), mode="RGB")
        photo = ImageTk.PhotoImage(pil)

        canvas.image = photo
        canvas.create_image(0, 0, anchor="nw", image=photo)

        def set_result(val):
            self.app.stop_work_panel()
            self._ok(val)

        def on_no():
            self.app.stop_work_panel()
            panel2 = self.app.start_work_panel(self.app._t("dlg.filter.problem_title"))

            tk.Label(panel2, text=self.app._t("dlg.filter.problem_ask"),
                     bg=self.app.BG_CARD, fg=self.app.FG,
                     font=("Segoe UI", 10, "bold"), justify="left", anchor="w", wraplength=500).pack(anchor="w", pady=(10, 20))

            def choose(val, warn_key):
                messagebox.showinfo(self.app._t("dlg.filter.problem_title"), self.app._t(warn_key))
                self.app.stop_work_panel()
                self._ok(val)

            if mode in ("cells", "soma"):
                RoundedButton(panel2, self.app._t("dlg.filter.prob_fused"),
                              command=lambda: choose("fused", "dlg.filter.warn_soma"),
                              width=500, height=40, radius=20,
                              bg=self.app.BG_SOFT, fg=self.app.FG,
                              hover_bg=self.app.WARN, hover_fg="#1a1c22",
                              active_bg=self.app.WARN).pack(pady=6, anchor="w")

            if mode in ("cells", "nuclei"):
                RoundedButton(panel2, self.app._t("dlg.filter.prob_artifacts"),
                              command=lambda: choose("artifacts", "dlg.filter.warn_slices" if getattr(self.app, "_deconv_was_on", False) else "dlg.filter.warn_slices_fallback"),
                              width=500, height=40, radius=20,
                              bg=self.app.BG_SOFT, fg=self.app.FG,
                              hover_bg=self.app.WARN, hover_fg="#1a1c22",
                              active_bg=self.app.WARN).pack(pady=6, anchor="w")

            if mode == "cells":
                RoundedButton(panel2, self.app._t("dlg.filter.prob_cut"),
                              command=lambda: choose("cut", "dlg.filter.warn_cut"),
                              width=500, height=40, radius=20,
                              bg=self.app.BG_SOFT, fg=self.app.FG,
                              hover_bg=self.app.WARN, hover_fg="#1a1c22",
                              active_bg=self.app.WARN).pack(pady=6, anchor="w")

            if mode in ("cells", "nuclei", "granules"):
                RoundedButton(panel2, self.app._t("dlg.filter.prob_not_cells"),
                              command=lambda: choose("delete", "dlg.filter.warn_delete"),
                              width=500, height=40, radius=20,
                              bg=self.app.BG_SOFT, fg=self.app.FG,
                              hover_bg=self.app.WARN, hover_fg="#1a1c22",
                              active_bg=self.app.WARN).pack(pady=6, anchor="w")

            self.app.work_panel_footer(
                panel2, on_ok=None, ok_text=None,
                cancel=True, cancel_value=None,
                cancel_text=self.app._t("btn.cancel"),
                help_text=""
            )
            for child in self.app._footer_btns.winfo_children():
                if isinstance(child, RoundedButton) and child._label == self.app._t("btn.cancel"):
                    child._command = lambda: set_result(None)

        self.app.work_panel_footer(
            panel, on_ok=lambda: set_result("ok"),
            ok_text=self.app._t("btn.yes"),
            extra=[(self.app._t("btn.no"), on_no)],
            help_text=""
        )
        for child in self.app._footer_btns.winfo_children():
            if isinstance(child, RoundedButton) and child._label == self.app._t("btn.no"):
                child._command = on_no

    def _dialog_cell_deletion(self, kw):
        mip, cell_mip, colors = kw["mip"], kw["cell_mip"], kw["colors"]
        y_dim, x_dim = mip.shape
        panel = self.app.start_work_panel(self.app._t("dlg.filter.delete_title"))

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

        def render():
            rgb = np.stack([base_small, base_small, base_small], axis=-1).astype(np.float32)
            for cid, small in masks_small.items():
                if not small.any():
                    continue
                if state[cid] == 1:
                    col = np.array([255, 60, 60], dtype=np.float32)
                    a = 0.8
                else:
                    col = np.array(colors.get(cid, (1.0, 0.5, 0.0)), dtype=np.float32) * 255.0
                    a = 0.45
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
            state[cid] = 1 if state[cid] == 0 else 0
            render()

        canvas.bind("<Button-1>", on_click)
        render()

        def finish():
            self.app.stop_work_panel()
            self._ok([cid for cid, st in state.items() if st == 1])

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            cancel=True, cancel_value=[],
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.filter.delete_help")
        )

    def _dialog_preprocess_choice(self, kw):
        need_ds_q = bool(kw.get("need_downsample_q", False))
        need_crop_q = bool(kw.get("need_crop_q", False))
        img_w = int(kw.get("img_w", 0))
        img_h = int(kw.get("img_h", 0))
        allow_ds_default = bool(kw.get("allow_downsample_default", True))
        do_crop_default = bool(kw.get("do_crop_default", False))

        panel = self.app.start_work_panel(self.app._t("dlg.preprocess.title"))

        tk.Label(
            panel,
            text=self.app._t("dlg.preprocess.size_line", w=img_w, h=img_h),
            bg=self.app.BG_CARD, fg=self.app.FG,
            font=("Segoe UI", 10, "bold"),
            anchor="w", justify="left"
        ).pack(anchor="w", pady=(0, 10))

        var_ds = tk.BooleanVar(value=allow_ds_default if need_ds_q else True)
        var_crop = tk.BooleanVar(value=do_crop_default if need_crop_q else False)

        if need_ds_q:
            chk_ds = tk.Checkbutton(
                panel,
                text=self.app._t("dlg.preprocess.ds_question"),
                variable=var_ds,
                bg=self.app.BG_CARD, fg=self.app.FG,
                activebackground=self.app.BG_CARD,
                activeforeground=self.app.ACCENT,
                selectcolor=self.app.BG_SOFT,
                relief="flat", bd=0,
                font=("Segoe UI", 10), cursor="hand2",
                highlightthickness=0,
                anchor="w", justify="left", wraplength=460,
            )
            chk_ds.pack(anchor="w", pady=(0, 8))

        if need_crop_q:
            chk_crop = tk.Checkbutton(
                panel,
                text=self.app._t("dlg.preprocess.crop_question"),
                variable=var_crop,
                bg=self.app.BG_CARD, fg=self.app.FG,
                activebackground=self.app.BG_CARD,
                activeforeground=self.app.ACCENT,
                selectcolor=self.app.BG_SOFT,
                relief="flat", bd=0,
                font=("Segoe UI", 10), cursor="hand2",
                highlightthickness=0,
                anchor="w", justify="left", wraplength=460,
            )
            chk_crop.pack(anchor="w", pady=(0, 8))

        def finish():
            self.app.stop_work_panel()
            self._ok({
                "allow_downsample": bool(var_ds.get()) if need_ds_q else True,
                "do_crop": bool(var_crop.get()) if need_crop_q else False,
            })

        def cancel():
            self.app.stop_work_panel()
            self._ok(None)

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            cancel=True, cancel_value=None,
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.preprocess.help"),
            on_cancel=cancel,
        )

    def _dialog_crop_editor(self, kw):
        mip = kw["mip"]
        xy_um = float(kw.get("xy_um", 1.0))
        if xy_um <= 0:
            xy_um = 1.0

        while mip.ndim > 2:
            mip = np.max(mip, axis=0)

        img_h, img_w = mip.shape[:2]

        panel = self.app.start_work_panel(self.app._t("dlg.crop.title"))

        min_size = 16
        min_size_um = min_size * xy_um

        state = {
            "x0": max(0, img_w // 8),
            "y0": max(0, img_h // 8),
            "x1": max(min_size, img_w - img_w // 8),
            "y1": max(min_size, img_h - img_h // 8),
            "drag": None,
            "drag_off_x": 0,
            "drag_off_y": 0,
        }

        if state["x1"] - state["x0"] < min_size:
            state["x0"] = 0
            state["x1"] = min(img_w, min_size)
        if state["y1"] - state["y0"] < min_size:
            state["y0"] = 0
            state["y1"] = min(img_h, min_size)

        info_row = tk.Frame(panel, bg=self.app.BG_CARD)
        info_row.pack(side="bottom", fill="x", pady=(8, 0))

        lbl_coords = tk.Label(
            info_row, text="",
            bg=self.app.BG_CARD, fg=self.app.FG_DIM,
            font=("Consolas", 9),
            anchor="w", justify="left"
        )
        lbl_coords.pack(side="left")

        size_row = tk.Frame(panel, bg=self.app.BG_CARD)
        size_row.pack(side="bottom", fill="x", pady=(6, 0))

        tk.Label(size_row, text=self.app._t("dlg.crop.width_um"),
                 bg=self.app.BG_CARD, fg=self.app.FG,
                 font=("Segoe UI", 9)).pack(side="left")

        var_w = tk.StringVar(value=f"{(state['x1'] - state['x0']) * xy_um:.2f}")
        ent_w = tk.Entry(size_row, textvariable=var_w, width=8,
                         bg=self.app.BG_SOFT, fg=self.app.FG,
                         insertbackground=self.app.ACCENT,
                         relief="flat", bd=0,
                         highlightthickness=1,
                         highlightbackground=self.app.BORDER,
                         highlightcolor=self.app.ACCENT,
                         font=("Consolas", 10))
        ent_w.pack(side="left", padx=(6, 16))

        tk.Label(size_row, text=self.app._t("dlg.crop.height_um"),
                 bg=self.app.BG_CARD, fg=self.app.FG,
                 font=("Segoe UI", 9)).pack(side="left")

        var_h = tk.StringVar(value=f"{(state['y1'] - state['y0']) * xy_um:.2f}")
        ent_h = tk.Entry(size_row, textvariable=var_h, width=8,
                         bg=self.app.BG_SOFT, fg=self.app.FG,
                         insertbackground=self.app.ACCENT,
                         relief="flat", bd=0,
                         highlightthickness=1,
                         highlightbackground=self.app.BORDER,
                         highlightcolor=self.app.ACCENT,
                         font=("Consolas", 10))
        ent_h.pack(side="left", padx=(6, 16))

        def apply_size_fields(*_):
            try:
                new_w_um = float(str(var_w.get()).replace(",", "."))
                new_h_um = float(str(var_h.get()).replace(",", "."))
            except Exception:
                return

            new_w = int(round(new_w_um / xy_um))
            new_h = int(round(new_h_um / xy_um))

            new_w = max(min_size, min(img_w, new_w))
            new_h = max(min_size, min(img_h, new_h))

            cx = (state["x0"] + state["x1"]) / 2.0
            cy = (state["y0"] + state["y1"]) / 2.0

            x0 = int(round(cx - new_w / 2.0))
            y0 = int(round(cy - new_h / 2.0))
            x0 = max(0, min(img_w - new_w, x0))
            y0 = max(0, min(img_h - new_h, y0))

            state["x0"] = x0
            state["x1"] = x0 + new_w
            state["y0"] = y0
            state["y1"] = y0 + new_h

            redraw()

        ent_w.bind("<Return>", apply_size_fields)
        ent_w.bind("<FocusOut>", apply_size_fields)
        ent_h.bind("<Return>", apply_size_fields)
        ent_h.bind("<FocusOut>", apply_size_fields)

        canvas = tk.Canvas(panel, bg="#101216", highlightthickness=0)
        canvas.pack(side="top", fill="both", expand=True, pady=(4, 8))

        view = {"s": 1.0, "x0": 0, "y0": 0, "photo": None}

        base_u8 = self._to_uint8(mip)
        base_pil = self._pil_from_gray(base_u8)

        def _handle_boxes():
            x0 = view["x0"] + state["x0"] * view["s"]
            x1 = view["x0"] + state["x1"] * view["s"]
            y0 = view["y0"] + state["y0"] * view["s"]
            y1 = view["y0"] + state["y1"] * view["s"]
            hs = 7
            handles = {
                "tl": (x0, y0), "tr": (x1, y0),
                "bl": (x0, y1), "br": (x1, y1),
                "t": ((x0 + x1) / 2, y0),
                "b": ((x0 + x1) / 2, y1),
                "l": (x0, (y0 + y1) / 2),
                "r": (x1, (y0 + y1) / 2),
            }
            return (x0, y0, x1, y1, hs, handles)

        def redraw(event=None):
            cw = canvas.winfo_width()
            ch = canvas.winfo_height()
            if cw <= 2 or ch <= 2:
                return

            s = min(cw / img_w, ch / img_h)
            new_w = max(1, int(img_w * s))
            new_h = max(1, int(img_h * s))
            ox = (cw - new_w) // 2
            oy = (ch - new_h) // 2

            view["s"] = s
            view["x0"] = ox
            view["y0"] = oy

            pil = base_pil.resize((new_w, new_h), Image.BILINEAR)
            photo = ImageTk.PhotoImage(pil)
            view["photo"] = photo

            canvas.delete("all")
            canvas.create_image(ox, oy, anchor="nw", image=photo)

            x0, y0, x1, y1, hs, handles = _handle_boxes()

            canvas.create_rectangle(x0, y0, x1, y1,
                                    outline="#ffd54a", width=2)
            for (_k, (hx, hy)) in handles.items():
                canvas.create_rectangle(hx - hs // 2, hy - hs // 2,
                                        hx + hs // 2, hy + hs // 2,
                                        outline="#ffd54a", width=2,
                                        fill="#1a1c22")

            lbl_coords.config(
                text=self.app._t(
                    "dlg.crop.coords",
                    x0=state["x0"] * xy_um, x1=state["x1"] * xy_um,
                    y0=state["y0"] * xy_um, y1=state["y1"] * xy_um,
                    w=(state["x1"] - state["x0"]) * xy_um,
                    h=(state["y1"] - state["y0"]) * xy_um,
                )
            )

        def _canvas_to_img(sx, sy):
            s = view["s"]
            ox = view["x0"]
            oy = view["y0"]
            if s <= 0:
                return None
            x = int(round((sx - ox) / s))
            y = int(round((sy - oy) / s))
            return x, y

        def _find_handle(sx, sy):
            x0, y0, x1, y1, hs, handles = _handle_boxes()
            for key, (hx, hy) in handles.items():
                if abs(sx - hx) <= hs and abs(sy - hy) <= hs:
                    return key
            if x0 <= sx <= x1 and y0 <= sy <= y1:
                return "move"
            return None

        def on_press(event):
            key = _find_handle(event.x, event.y)
            if key is None:
                return
            state["drag"] = key
            if key == "move":
                pt = _canvas_to_img(event.x, event.y)
                if pt is not None:
                    state["drag_off_x"] = pt[0] - state["x0"]
                    state["drag_off_y"] = pt[1] - state["y0"]

        def on_motion(event):
            key = _find_handle(event.x, event.y)
            cursors = {
                "tl": "top_left_corner", "tr": "top_right_corner",
                "bl": "bottom_left_corner", "br": "bottom_right_corner",
                "t": "top_side", "b": "bottom_side",
                "l": "left_side", "r": "right_side",
                "move": "fleur",
            }
            if key is not None:
                canvas.configure(cursor=cursors.get(key, "hand2"))
            else:
                canvas.configure(cursor="")

            if state["drag"] is None:
                return

            pt = _canvas_to_img(event.x, event.y)
            if pt is None:
                return
            mx, my = pt

            drag = state["drag"]
            x0, y0, x1, y1 = state["x0"], state["y0"], state["x1"], state["y1"]

            if drag == "move":
                new_w = x1 - x0
                new_h = y1 - y0
                nx0 = mx - state["drag_off_x"]
                ny0 = my - state["drag_off_y"]
                nx0 = max(0, min(img_w - new_w, nx0))
                ny0 = max(0, min(img_h - new_h, ny0))
                state["x0"] = nx0
                state["x1"] = nx0 + new_w
                state["y0"] = ny0
                state["y1"] = ny0 + new_h
            else:
                if "l" in drag:
                    x0 = max(0, min(x1 - min_size, mx))
                if "r" in drag:
                    x1 = min(img_w, max(x0 + min_size, mx))
                if "t" in drag:
                    y0 = max(0, min(y1 - min_size, my))
                if "b" in drag:
                    y1 = min(img_h, max(y0 + min_size, my))
                state["x0"], state["x1"] = int(x0), int(x1)
                state["y0"], state["y1"] = int(y0), int(y1)

            var_w.set(f"{(state['x1'] - state['x0']) * xy_um:.2f}")
            var_h.set(f"{(state['y1'] - state['y0']) * xy_um:.2f}")
            redraw()

        def on_release(event):
            state["drag"] = None

        canvas.bind("<Configure>", redraw)
        canvas.bind("<ButtonPress-1>", on_press)
        canvas.bind("<B1-Motion>", on_motion)
        canvas.bind("<ButtonRelease-1>", on_release)
        canvas.bind("<Motion>", lambda e: on_motion(e) if state["drag"] is None else None)

        def finish():
            self.app.stop_work_panel()
            self._ok((state["y0"], state["y1"], state["x0"], state["x1"]))

        def cancel():
            self.app.stop_work_panel()
            self._ok(None)

        canvas.after(30, redraw)

        self.app.work_panel_footer(
            panel, on_ok=finish,
            ok_text=self.app._t("btn.apply"),
            cancel=True, cancel_value=None,
            cancel_text=self.app._t("btn.cancel"),
            help_text=self.app._t("dlg.crop.help"),
            on_cancel=cancel,
        )

class RoleDropdown(tk.Frame):
    ORDER = ("cells", "dapi", "granules", "none")

    def __init__(self, parent, app, on_change,
                 labels, taken_label,
                 width=180, height=30):
        super().__init__(parent, bg=app.BG_CARD)

        self._app = app
        self._on_change = on_change
        self._labels = labels
        self._taken_label = taken_label
        self._state = "none"
        self._used = {"cells": None, "dapi": None, "granules": None}
        self._index = 0
        self._btn_w = width
        self._btn_h = height

        self._btn = tk.Canvas(self, width=self._btn_w, height=self._btn_h,
                              bg=app.BG_CARD, highlightthickness=0, bd=0)
        self._btn.pack(side="top")
        self._btn.bind("<Button-1>", self._on_click)
        self._btn.bind("<Enter>", lambda e: self._draw(hover=True))
        self._btn.bind("<Leave>", lambda e: self._draw(hover=False))
        self._hover = False
        self._draw()

    def _color_for(self, role):
        return {
            "cells": self._app.OK,
            "dapi": self._app.ACCENT,
            "granules": self._app.WARN,
            "none": self._app.FG_MUTED,
        }.get(role, self._app.FG)

    def _rounded(self, x1, y1, x2, y2, r, **kw):
        pts = [
            x1 + r, y1, x2 - r, y1, x2, y1,
            x2, y1 + r, x2, y2 - r, x2, y2,
            x2 - r, y2, x1 + r, y2, x1, y2,
            x1, y2 - r, x1, y1 + r, x1, y1,
            ]
        return self._btn.create_polygon(pts, smooth=True, **kw)

    def _draw(self, hover=None):
        if hover is not None:
            self._hover = hover
        self._btn.delete("all")
        bg = self._app.BG_HOVER if self._hover else self._app.BG_SOFT
        fg = self._color_for(self._state)
        r = self._btn_h // 2
        self._rounded(1, 1, self._btn_w - 1, self._btn_h - 1, r,
                      fill=bg, outline=self._app.BORDER, width=1)
        self._btn.create_text(12, self._btn_h // 2, anchor="w",
                              text=self._labels.get(self._state, "—"),
                              fill=fg, font=("Segoe UI", 9, "bold"))
        self._btn.create_text(self._btn_w - 14, self._btn_h // 2, anchor="e",
                              text="▾", fill=self._app.FG,
                              font=("Segoe UI", 11, "bold"))

    def set_state(self, role, used, index):
        self._state = role
        self._used = used
        self._index = index
        self._draw()

    def _on_click(self, event):
        menu = tk.Menu(self, tearoff=0,
                       bg=self._app.BG_SOFT, fg=self._app.FG,
                       activebackground=self._app.ACCENT,
                       activeforeground="#1a1c22",
                       relief="flat", bd=0,
                       font=("Segoe UI", 10))

        for role in self.ORDER:
            label_base = self._labels.get(role, role)
            taken_by = self._used.get(role) if role != "none" else None
            is_taken_elsewhere = (taken_by is not None and taken_by != self._index)
            is_current = (self._state == role)
            prefix = "✓  " if is_current else "     "
            if is_taken_elsewhere:
                menu.add_command(
                    label=f"{prefix}{label_base}   —   ({self._taken_label})",
                    state="disabled",
                )
            else:
                menu.add_command(
                    label=f"{prefix}{label_base}",
                    command=lambda r=role: self._choose(r),
                )

        try:
            menu.tk_popup(event.x_root, event.y_root)
        finally:
            menu.grab_release()

    def _choose(self, role):
        if role == self._state:
            return
        self._on_change(role)

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
        self._enabled = True

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
        if not self._enabled:
            bg = self._bg
            fg = "#5a5e68"
        else:
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
        if not self._enabled:
            return
        self._hover = True
        self.configure(cursor="hand2")
        self._draw()

    def _on_leave(self, _):
        self._hover = False
        self.configure(cursor="")
        self._draw()

    def _on_click(self, _):
        if not self._enabled:
            return
        if self._var is not None:
            self._var.set(self._value)

    def set_label(self, text):
        self._label = text
        self._measure()
        self._draw()

    def set_enabled(self, enabled: bool):
        self._enabled = bool(enabled)
        if not self._enabled:
            self._hover = False
            self.configure(cursor="")
        else:
            self.configure(cursor="hand2")
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
        self._pending_metrics_config = default_metrics_config()
        self._pending_channels_config = None
        self._current_channel_roles = {}
        self._active_channel_finish = None
        self._xy_from_metadata = False
        self._z_from_metadata = False

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
        self._has_started = False

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

        for lbl, role in getattr(self, "_src_labels", []):
            try:
                if role == "cells":
                    if hasattr(self, "_last_optics_state") and "cells" in self._last_optics_state:
                        self._apply_optics_label(
                            lbl, self._last_optics_state["cells"],
                            entry_na=self.entry_na_cells, entry_wl=self.entry_wl_cells,
                            var_na=self.var_na_cells, var_wl=self.var_wl_cells
                        )
                    else:
                        lbl.config(text=self._t("ui.settings.default"), fg=self.FG_DIM)
                elif role == "dapi":
                    if hasattr(self, "_last_optics_state") and "dapi" in self._last_optics_state:
                        self._apply_optics_label(
                            lbl, self._last_optics_state["dapi"],
                            entry_na=self.entry_na_dapi, entry_wl=self.entry_wl_dapi,
                            var_na=self.var_na_dapi, var_wl=self.var_wl_dapi
                        )
                    else:
                        lbl.config(text=self._t("ui.settings.default"), fg=self.FG_DIM)
                elif role == "granules":
                    if hasattr(self, "_last_optics_state") and "granules" in self._last_optics_state:
                        self._apply_optics_label(
                            lbl, self._last_optics_state["granules"],
                            entry_na=self.entry_na_granules, entry_wl=self.entry_wl_granules,
                            var_na=self.var_na_granules, var_wl=self.var_wl_granules
                        )
                    else:
                        lbl.config(text=self._t("ui.settings.default"), fg=self.FG_DIM)
                elif role == "xy":
                    if getattr(self, "_xy_from_metadata", False):
                        lbl.config(text=self._t("ui.settings.from_metadata"), fg=self.OK)
                    else:
                        lbl.config(text=self._t("ui.settings.default"), fg=self.FG_DIM)
                elif role == "z":
                    if getattr(self, "_z_from_metadata", False):
                        lbl.config(text=self._t("ui.settings.from_metadata"), fg=self.OK)
                    else:
                        lbl.config(text=self._t("ui.settings.default"), fg=self.FG_DIM)
                else:
                    lbl.config(text=self._t(role))
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

    def _build_welcome(self, parent):
        self.welcome_container = tk.Frame(parent, bg=self.BG)
        self.welcome_container.pack(fill="both", expand=True)

        wrap = tk.Frame(self.welcome_container, bg=self.BG_CARD,
                        bd=0, highlightthickness=1,
                        highlightbackground=self.BORDER,
                        highlightcolor=self.BORDER)
        wrap.pack(fill="both", expand=True)

        header = tk.Frame(wrap, bg=self.BG_CARD)
        header.pack(fill="x", padx=16, pady=(12, 2))

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

        self.lbl_main = tk.Label(header,
                                 text=self._t("ui.settings.title_main"),
                                 bg=self.BG_CARD, fg=self.FG,
                                 font=("Segoe UI", 10, "bold"),
                                 anchor="w")
        self.lbl_main.pack(fill="x")

        self._settings_static_labels.append((self.lbl_main, "ui.settings.title_main"))

        body = tk.Frame(wrap, bg=self.BG_CARD)
        body.pack(fill="x", padx=16, pady=(6, 12))

        self._build_settings(body)

    def _switch_to_welcome_mode(self):
        if not self._analysis_mode:
            return
        self._analysis_mode = False
        self.btn_back_welcome.pack_forget()
        self.log_container.pack_forget()
        self.welcome_container.pack(fill="both", expand=True)
        self.btn_to_log.pack(side="right")

        if self._has_started:
            self.lbl_main.pack_forget()

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

        row = self._row(s, "ui.settings.auto_soma")
        self.var_auto_soma = tk.BooleanVar(value=False)
        chk_auto_soma = self._check(row, self._t("ui.settings.auto_soma_check"),
                                    self.var_auto_soma)
        chk_auto_soma.pack(side="left")
        self._settings_checks.append((chk_auto_soma, "ui.settings.auto_soma_check"))

        row = self._row(s, "ui.settings.xy")
        self.var_xy = tk.StringVar(value="0.144")
        self.entry_xy = self._entry(row, self.var_xy, width=10)
        self.entry_xy.pack(side="left")
        self.lbl_xy_src = tk.Label(row, text=self._t("ui.settings.default"),
                                   bg=self.BG_CARD, fg=self.FG_DIM,
                                   font=("Segoe UI", 9, "italic"))
        self.lbl_xy_src.pack(side="left", padx=(8, 0))
        self._src_labels.append((self.lbl_xy_src, "xy"))

        row = self._row(s, "ui.settings.z")
        self.var_z = tk.StringVar(value="1.0")
        self.entry_z = self._entry(row, self.var_z, width=10)
        self.entry_z.pack(side="left")
        self.lbl_z_src = tk.Label(row, text=self._t("ui.settings.default"),
                                  bg=self.BG_CARD, fg=self.FG_DIM,
                                  font=("Segoe UI", 9, "italic"))
        self.lbl_z_src.pack(side="left", padx=(8, 0))
        self._src_labels.append((self.lbl_z_src, "z"))

        hint = tk.Label(s, text=self._t("ui.settings.scale_hint"),
                        bg=self.BG_CARD, fg=self.FG_DIM,
                        font=("Segoe UI", 9, "italic"),
                        anchor="w", justify="left",
                        wraplength=760)
        hint.pack(fill="x", pady=(0, 6))
        self._settings_static_labels.append((hint, "ui.settings.scale_hint"))

        row = self._row(s, "ui.settings.cells", width=14)
        self.lbl_na_cells = tk.Label(row, text=self._t("ui.settings.na"),
                                     bg=self.BG_CARD, fg=self.FG,
                                     font=("Segoe UI", 10))
        self.lbl_na_cells.pack(side="left", padx=(0, 4))
        self._settings_labels.append((self.lbl_na_cells, "ui.settings.na"))
        self.var_na_cells = tk.StringVar(value="0.95")
        self.entry_na_cells = self._entry(row, self.var_na_cells, width=6)
        self.entry_na_cells.pack(side="left")
        self.lbl_wl_cells = tk.Label(row, text=self._t("ui.settings.lambda_nm"),
                                     bg=self.BG_CARD, fg=self.FG,
                                     font=("Segoe UI", 10))
        self.lbl_wl_cells.pack(side="left", padx=(14, 4))
        self._settings_labels.append((self.lbl_wl_cells, "ui.settings.lambda_nm"))
        self.var_wl_cells = tk.StringVar(value="699.5")
        self.entry_wl_cells = self._entry(row, self.var_wl_cells, width=8)
        self.entry_wl_cells.pack(side="left")
        self.lbl_eps_cells = tk.Label(row, text=self._t("ui.settings.epsilon"),
                                      bg=self.BG_CARD, fg=self.FG,
                                      font=("Segoe UI", 10))
        self.lbl_eps_cells.pack(side="left", padx=(14, 4))
        self._settings_labels.append((self.lbl_eps_cells, "ui.settings.epsilon"))
        self.var_eps_cells = tk.StringVar(value="0.005")
        self.entry_eps_cells = self._entry(row, self.var_eps_cells, width=8)
        self.entry_eps_cells.pack(side="left")
        self.lbl_cells_src = tk.Label(row, text=self._t("ui.settings.default"),
                                      bg=self.BG_CARD, fg=self.FG_DIM,
                                      font=("Segoe UI", 9, "italic"))
        self.lbl_cells_src.pack(side="left", padx=(12, 0))
        self._src_labels.append((self.lbl_cells_src, "cells"))

        row = self._row(s, "ui.settings.dapi", width=14)
        self.lbl_na_dapi = tk.Label(row, text=self._t("ui.settings.na"),
                                    bg=self.BG_CARD, fg=self.FG,
                                    font=("Segoe UI", 10))
        self.lbl_na_dapi.pack(side="left", padx=(0, 4))
        self._settings_labels.append((self.lbl_na_dapi, "ui.settings.na"))
        self.var_na_dapi = tk.StringVar(value="0.95")
        self.entry_na_dapi = self._entry(row, self.var_na_dapi, width=6)
        self.entry_na_dapi.pack(side="left")
        self.lbl_wl_dapi = tk.Label(row, text=self._t("ui.settings.lambda_nm"),
                                    bg=self.BG_CARD, fg=self.FG,
                                    font=("Segoe UI", 10))
        self.lbl_wl_dapi.pack(side="left", padx=(14, 4))
        self._settings_labels.append((self.lbl_wl_dapi, "ui.settings.lambda_nm"))
        self.var_wl_dapi = tk.StringVar(value="460")
        self.entry_wl_dapi = self._entry(row, self.var_wl_dapi, width=8)
        self.entry_wl_dapi.pack(side="left")
        self.lbl_eps_dapi = tk.Label(row, text=self._t("ui.settings.epsilon"),
                                     bg=self.BG_CARD, fg=self.FG,
                                     font=("Segoe UI", 10))
        self.lbl_eps_dapi.pack(side="left", padx=(14, 4))
        self._settings_labels.append((self.lbl_eps_dapi, "ui.settings.epsilon"))
        self.var_eps_dapi = tk.StringVar(value="0.005")
        self.entry_eps_dapi = self._entry(row, self.var_eps_dapi, width=8)
        self.entry_eps_dapi.pack(side="left")
        self.lbl_dapi_src = tk.Label(row, text=self._t("ui.settings.default"),
                                     bg=self.BG_CARD, fg=self.FG_DIM,
                                     font=("Segoe UI", 9, "italic"))
        self.lbl_dapi_src.pack(side="left", padx=(12, 0))
        self._src_labels.append((self.lbl_dapi_src, "dapi"))

        row = self._row(s, "ui.settings.granules", width=14)
        self.lbl_na_granules = tk.Label(row, text=self._t("ui.settings.na"),
                                        bg=self.BG_CARD, fg=self.FG,
                                        font=("Segoe UI", 10))
        self.lbl_na_granules.pack(side="left", padx=(0, 4))
        self._settings_labels.append((self.lbl_na_granules, "ui.settings.na"))
        self.var_na_granules = tk.StringVar(value="0.95")
        self.entry_na_granules = self._entry(row, self.var_na_granules, width=6)
        self.entry_na_granules.pack(side="left")
        self.lbl_wl_granules = tk.Label(row, text=self._t("ui.settings.lambda_nm"),
                                        bg=self.BG_CARD, fg=self.FG,
                                        font=("Segoe UI", 10))
        self.lbl_wl_granules.pack(side="left", padx=(14, 4))
        self._settings_labels.append((self.lbl_wl_granules, "ui.settings.lambda_nm"))
        self.var_wl_granules = tk.StringVar(value="520")
        self.entry_wl_granules = self._entry(row, self.var_wl_granules, width=8)
        self.entry_wl_granules.pack(side="left")
        self.lbl_eps_granules = tk.Label(row, text=self._t("ui.settings.epsilon"),
                                         bg=self.BG_CARD, fg=self.FG,
                                         font=("Segoe UI", 10))
        self.lbl_eps_granules.pack(side="left", padx=(14, 4))
        self._settings_labels.append((self.lbl_eps_granules, "ui.settings.epsilon"))
        self.var_eps_granules = tk.StringVar(value="0.005")
        self.entry_eps_granules = self._entry(row, self.var_eps_granules, width=8)
        self.entry_eps_granules.pack(side="left")
        self.lbl_granules_src = tk.Label(row, text=self._t("ui.settings.default"),
                                         bg=self.BG_CARD, fg=self.FG_DIM,
                                         font=("Segoe UI", 9, "italic"))
        self.lbl_granules_src.pack(side="left", padx=(12, 0))
        self._src_labels.append((self.lbl_granules_src, "granules"))

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
            if isinstance(c, (tk.Checkbutton, tk.Entry, tk.Radiobutton,
                              Chip, RoundedButton)):
                self._settings_widgets.append(c)
            self._collect_settings_widgets(c)

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
                          extra=None, help_text=None, cancel_text=None,
                          on_cancel=None):
        if help_text:
            self._footer_help.config(text=help_text)

        if cancel_text is None:
            cancel_text = self._t("btn.cancel")

        def _cancel():
            self.stop_work_panel()
            if on_cancel is not None:
                on_cancel()
            else:
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
                self._pending_metrics_config = default_metrics_config()
                self._pending_channels_config = None
        else:
            d = filedialog.askdirectory()
            if d:
                self.var_file_path.set(d)
                self._pending_metrics_config = default_metrics_config()
                self._pending_channels_config = None

        path = self.var_file_path.get().strip()
        if not path or not os.path.exists(path):
            return

        if hasattr(self, "btn_browse"):
            self.btn_browse.config_state("disabled")
            self.btn_browse.set_text("Загрузка...")

        import threading
        threading.Thread(target=self._preload_metadata_thread, args=(path,), daemon=True).start()


    def _apply_optics_label(self, lbl_src, ch, entry_na=None, entry_wl=None,
                            var_na=None, var_wl=None):
        if ch is None:
            if entry_na is not None:
                entry_na.config(fg=self.FG_MUTED)
            if entry_wl is not None:
                entry_wl.config(fg=self.FG_MUTED)
            lbl_src.config(text=self._t("ui.settings.default"), fg=self.FG_DIM)
            return False

        missing = []
        found = False

        if ch.get('na') is not None:
            if var_na is not None:
                var_na.set(f"{ch['na']:.3f}")
            if entry_na is not None:
                entry_na.config(fg=self.FG)
            found = True
        else:
            if entry_na is not None:
                entry_na.config(fg=self.WARN)
            missing.append(self._t("ui.settings.meta_na"))

        if ch.get('emission_um') is not None:
            if var_wl is not None:
                var_wl.set(f"{ch['emission_um']*1000:.1f}")
            if entry_wl is not None:
                entry_wl.config(fg=self.FG)
            found = True
        else:
            if entry_wl is not None:
                entry_wl.config(fg=self.WARN)
            missing.append(self._t("ui.settings.meta_wl"))

        if missing:
            missing_str = ", ".join(missing)
            lbl_src.config(
                text=self._t("ui.settings.not_in_metadata", missing=missing_str),
                fg=self.WARN
            )
        elif found:
            lbl_src.config(text=self._t("ui.settings.from_metadata"), fg=self.OK)
        else:
            lbl_src.config(text=self._t("ui.settings.default"), fg=self.FG_DIM)

        return found

    def _preload_metadata_thread(self, path):
        try:
            sample = path
            if os.path.isdir(path):
                exts = ['.tif', '.tiff', '.nd2', '.ome.tif', '.ome.tiff', '.oib', '.oif', '.oir']
                files = sorted(
                    os.path.join(path, f) for f in os.listdir(path)
                    if any(f.lower().endswith(e) for e in exts))
                if not files:
                    return
                sample = files[0]

            self._current_filepath = sample

            ctx = AnalysisContext()
            ctx.bridge = self.bridge
            ctx.lang = self.var_lang.get()
            file_io = FileIO(ctx)
            optics = Optics(ctx, file_io)

            meta = file_io._extract_voxel_size_from_file(sample)

            try:
                op = optics._extract_optics_from_file(sample)
                channels = op.get('channels', [])
            except Exception:
                channels = []

            def _update_base_meta():
                if meta.get('xy_um'):
                    self.var_xy.set(f"{meta['xy_um']:.4f}")
                    self.lbl_xy_src.config(text=self._t("ui.settings.from_metadata"), fg=self.OK)
                    self._xy_from_metadata = True
                else:
                    self._xy_from_metadata = False

                if meta.get('z_um'):
                    self.var_z.set(f"{meta['z_um']:.4f}")
                    self.lbl_z_src.config(text=self._t("ui.settings.from_metadata"), fg=self.OK)
                    self._z_from_metadata = True
                else:
                    self._z_from_metadata = False

                if channels:
                    cell_idx = op.get('astro_channel_index', 0)
                    dapi_idx = op.get('dapi_channel_index')
                    self._update_optics_for_channels(sample, cell_idx, dapi_idx, None)

            self.root.after(0, _update_base_meta)

            img = file_io._read_image_universal(sample)
            if img is not None:
                if img.ndim == 2:
                    img = img[np.newaxis, np.newaxis, ...]
                    self._probe_is_2d = True
                elif img.ndim == 3:
                    if img.shape[0] <= 8 and img.shape[1] > img.shape[0]:
                        img = img[:, np.newaxis, ...]
                        self._probe_is_2d = True
                    else:
                        img = img[np.newaxis, ...]
                        self._probe_is_2d = (img.shape[1] == 1)
                elif img.ndim == 4:
                    if img.shape[0] <= 8:
                        pass
                    else:
                        img = np.transpose(img, (1, 0, 2, 3))
                    self._probe_is_2d = (img.shape[1] == 1)

                def norm(m):
                    v = m[m > 0]
                    if len(v) == 0:
                        return m
                    p0, p99 = np.percentile(v, [1, 99])
                    return np.clip((m - p0) / (p99 - p0 + 1e-6), 0, 1)

                mips = [norm(np.max(img[i], axis=0)) for i in range(img.shape[0])]
                auto_cell_idx = op.get('astro_channel_index', 0) if 'op' in locals() else 0
                auto_dapi_idx = op.get('dapi_channel_index') if 'op' in locals() else None

                res = self.bridge._ask(
                    "channel_selector",
                    mips=mips,
                    cell_idx_default=auto_cell_idx,
                    dapi_idx_default=auto_dapi_idx,
                    has_dapi_default=(auto_dapi_idx is not None),
                )

                if res:
                    self._pending_channels_config = res
                    self._pending_metrics_config = res.get("metrics_config",
                                                           default_metrics_config(is_2d=getattr(self, "_probe_is_2d", False)))
                    self.root.after(0, lambda: self._update_optics_for_channels(
                        sample, res.get("cell_idx"), res.get("dapi_idx"), res.get("granule_idx")
                    ))

        except Exception as e:
            self.append_log(self._t("log.error_prefix", err=e))
            self.append_log(traceback.format_exc())
        finally:
            self.root.after(0, lambda: self._reset_browse_btn())

    def _reset_browse_btn(self):
        if hasattr(self, "btn_browse"):
            self.btn_browse.config_state("normal")
            self.btn_browse.set_text(self._t("ui.settings.browse"))

    def _update_optics_for_channels(self, sample_path, cell_idx, dapi_idx, granule_idx):
        if not sample_path or not os.path.exists(sample_path):
            return

        if hasattr(self, 'analyzer'):
            ctx = self.analyzer.ctx
        else:
            ctx = AnalysisContext()
            ctx.bridge = self.bridge
            ctx.lang = self.var_lang.get()

        file_io = FileIO(ctx)
        optics = Optics(ctx, file_io)

        try:
            op = optics._extract_optics_from_file(sample_path)
            channels = op.get('channels', [])
        except Exception as e:
            self.append_log(f"Optics error: {e}")
            return

        if not channels:
            return

        def _get_ch(idx):
            if idx is not None and 0 <= idx < len(channels):
                return channels[idx]
            return None

        ch_cells = _get_ch(cell_idx)
        ch_dapi = _get_ch(dapi_idx)
        ch_granules = _get_ch(granule_idx)

        if not hasattr(self, "_last_optics_state"):
            self._last_optics_state = {}

        if ch_cells is not None or ch_dapi is not None or ch_granules is not None:
            self._last_optics_state = {
                "cells": ch_cells,
                "dapi": ch_dapi,
                "granules": ch_granules
            }
        else:
            ch_cells = self._last_optics_state.get("cells")
            ch_dapi = self._last_optics_state.get("dapi")
            ch_granules = self._last_optics_state.get("granules")

        if ch_cells is not None:
            self._apply_optics_label(
                self.lbl_cells_src, ch_cells,
                entry_na=self.entry_na_cells, entry_wl=self.entry_wl_cells,
                var_na=self.var_na_cells, var_wl=self.var_wl_cells,
            )
        if ch_dapi is not None:
            self._apply_optics_label(
                self.lbl_dapi_src, ch_dapi,
                entry_na=self.entry_na_dapi, entry_wl=self.entry_wl_dapi,
                var_na=self.var_na_dapi, var_wl=self.var_wl_dapi,
            )
        if ch_granules is not None:
            self._apply_optics_label(
                self.lbl_granules_src, ch_granules,
                entry_na=self.entry_na_granules, entry_wl=self.entry_wl_granules,
                var_na=self.var_na_granules, var_wl=self.var_wl_granules,
            )

        self._update_settings_state_by_channels(dapi_idx is not None, granule_idx is not None)
        self.root.update_idletasks()

    def _read_metadata(self):
        self.append_log(self._t("log.metadata_reading", name="..."))
        self.root.update_idletasks()

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
                self._xy_from_metadata = True
                self.append_log(self._t("log.metadata_xy_ok", v=meta['xy_um']))
            else:
                self.lbl_xy_src.config(text=self._t("ui.settings.default"),
                                       fg=self.FG_DIM)
                self._xy_from_metadata = False
                self.append_log(self._t("log.metadata_xy_none"))

            if meta.get('z_um'):
                self.var_z.set(f"{meta['z_um']:.4f}")
                self.lbl_z_src.config(text=self._t("ui.settings.from_metadata"),
                                      fg=self.OK)
                self._z_from_metadata = True
                self.append_log(self._t("log.metadata_z_ok", v=meta['z_um']))
            else:
                self.lbl_z_src.config(text=self._t("ui.settings.default"),
                                      fg=self.FG_DIM)
                self._z_from_metadata = False
                self.append_log(self._t("log.metadata_z_none"))

            try:
                op = optics._extract_optics_from_file(sample)
                channels = op.get('channels', [])
            except Exception as e:
                self.append_log(self._t("log.metadata_optics_read_error", err=e))
                channels = []

            if channels:
                cell_idx = None
                dapi_idx = None
                granule_idx = None
                has_dapi = False
                has_granules = False

                roles_to_use = getattr(self, "_current_channel_roles", {})
                if not roles_to_use and getattr(self, "_pending_channels_config", None):
                    cfg = self._pending_channels_config
                    if cfg.get("cell_idx") is not None:
                        roles_to_use[cfg["cell_idx"]] = "cells"
                    if cfg.get("dapi_idx") is not None:
                        roles_to_use[cfg["dapi_idx"]] = "dapi"
                    if cfg.get("granule_idx") is not None:
                        roles_to_use[cfg["granule_idx"]] = "granules"

                if roles_to_use:
                    for i, r in roles_to_use.items():
                        if r == "cells": cell_idx = i
                        elif r == "dapi":
                            dapi_idx = i
                            has_dapi = True
                        elif r == "granules":
                            granule_idx = i
                            has_granules = True
                else:
                    cell_idx = op.get('astro_channel_index', 0)
                    dapi_idx = op.get('dapi_channel_index')
                    has_dapi = (dapi_idx is not None)

                if getattr(self, "_pending_channels_config", None):
                    self._pending_channels_config.update({
                        "cell_idx": cell_idx if cell_idx is not None else 0,
                        "has_dapi": has_dapi,
                        "dapi_idx": dapi_idx,
                        "has_granules": has_granules,
                        "granule_idx": granule_idx,
                    })
                else:
                    self._pending_channels_config = {
                        "cell_idx": cell_idx if cell_idx is not None else 0,
                        "has_dapi": has_dapi,
                        "dapi_idx": dapi_idx,
                        "has_granules": has_granules,
                        "granule_idx": granule_idx,
                        "metrics_config": dict(getattr(self, "_pending_metrics_config", {}))
                    }

                self._update_settings_state_by_channels(has_dapi, has_granules)
                self._update_optics_for_channels(sample, cell_idx, dapi_idx, granule_idx)

            self.append_log(self._t("log.metadata_fields_updated"))
            self.root.update_idletasks()
            StyledDialog.show(self.root,
                              self._t("ui.settings.metadata"),
                              self._t("log.metadata_fields_updated"),
                              kind="info")

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

    def _set_settings_state(self, state):
        enabled = (state != "disabled")
        for c in self._settings_widgets:
            try:
                if isinstance(c, RoundedButton):
                    c.config_state("normal" if enabled else "disabled")
                elif isinstance(c, Chip):
                    c.set_enabled(enabled)
                else:
                    c.configure(state=state)
            except Exception:
                pass

        for btn, _key in getattr(self, "_settings_buttons", []):
            if btn is getattr(self, "btn_run", None):
                continue
            if btn is getattr(self, "btn_exit", None):
                continue
            try:
                btn.config_state("normal" if enabled else "disabled")
            except Exception:
                pass

    def _update_settings_state_by_channels(self, has_dapi, has_granules):
        if not hasattr(self, "entry_na_dapi") or not hasattr(self, "entry_na_granules"):
            return

        is_locked = getattr(self, "_has_started", False)

        def _set_entry_state(widget, enabled):
            if enabled and not is_locked:
                widget.config(state="normal", disabledbackground=self.BG_SOFT, disabledforeground=self.FG)
            else:
                widget.config(state="disabled", disabledbackground=self.BG_CARD, disabledforeground=self.FG_MUTED)

        _set_entry_state(self.entry_na_dapi, has_dapi)
        _set_entry_state(self.entry_wl_dapi, has_dapi)
        _set_entry_state(self.entry_eps_dapi, has_dapi)

        if has_dapi and not is_locked:
            for w in (self.lbl_na_dapi, self.lbl_wl_dapi, self.lbl_eps_dapi):
                w.config(fg=self.FG)
            self.lbl_dapi_src.config(fg=self.OK)
        else:
            for w in (self.lbl_na_dapi, self.lbl_wl_dapi, self.lbl_eps_dapi, self.lbl_dapi_src):
                w.config(fg=self.FG_MUTED)

        _set_entry_state(self.entry_na_granules, has_granules)
        _set_entry_state(self.entry_wl_granules, has_granules)
        _set_entry_state(self.entry_eps_granules, has_granules)

        if has_granules and not is_locked:
            for w in (self.lbl_na_granules, self.lbl_wl_granules, self.lbl_eps_granules):
                w.config(fg=self.FG)
            self.lbl_granules_src.config(fg=self.OK)
        else:
            for w in (self.lbl_na_granules, self.lbl_wl_granules, self.lbl_eps_granules, self.lbl_granules_src):
                w.config(fg=self.FG_MUTED)

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

        self._has_started = True

        self._set_settings_state("disabled")

        ch_cfg = getattr(self, "_pending_channels_config", None)
        has_dapi = ch_cfg.get("has_dapi", False) if ch_cfg else False
        has_granules = ch_cfg.get("has_granules", False) if ch_cfg else False
        self._update_settings_state_by_channels(has_dapi, has_granules)

        if self._work_panel_active:
            if getattr(self, "_active_channel_finish", None) is not None:
                self._active_channel_finish()
                self._active_channel_finish = None
            else:
                self.stop_work_panel()
                self.bridge._ok(None)

        self._switch_to_analysis_mode()
        self.btn_run.config_state("disabled")

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
            "na_granules": f(self.var_na_granules.get()),
            "wl_granules": f(self.var_wl_granules.get()),
            "eps_granules": f(self.var_eps_granules.get()),
            "do_filter": bool(self.var_do_filter.get()),
            "auto_soma": bool(self.var_auto_soma.get()),
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

            self._current_filepath = files[0]

            analyzer = Analyzer()
            analyzer.ctx.bridge = self.bridge
            analyzer.ctx.lang = self.var_lang.get()
            analyzer.ctx.app = self

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

            analyzer.ctx._deconv_params['na'] = cfg["na_cells"]
            analyzer.ctx._deconv_params['emission_wavelength_um'] = \
                cfg["wl_cells"] / 1000.0
            analyzer.ctx._deconv_params['wiener_epsilon'] = cfg["eps_cells"]

            analyzer.ctx._dapi_deconv_params['na'] = cfg["na_dapi"]
            analyzer.ctx._dapi_deconv_params['emission_wavelength_um'] = \
                cfg["wl_dapi"] / 1000.0
            analyzer.ctx._dapi_deconv_params['wiener_epsilon'] = cfg["eps_dapi"]

            analyzer.ctx._granule_deconv_params['na'] = cfg["na_granules"]
            analyzer.ctx._granule_deconv_params['emission_wavelength_um'] = \
                cfg["wl_granules"] / 1000.0
            analyzer.ctx._granule_deconv_params['wiener_epsilon'] = cfg["eps_granules"]

            analyzer.ctx._gui_do_filter = cfg["do_filter"]
            analyzer.ctx.auto_soma = cfg["auto_soma"]

            mc = getattr(self, "_pending_metrics_config", None)
            if isinstance(mc, dict) and mc:
                analyzer.ctx.metrics_config = dict(mc)
            else:
                analyzer.ctx.metrics_config = default_metrics_config(
                    is_2d=getattr(self, "_probe_is_2d", False))

            ch_cfg = getattr(self, "_pending_channels_config", None)
            if ch_cfg and "cell_idx" in ch_cfg:
                analyzer.ctx.cell_channel_idx = ch_cfg.get("cell_idx", 0)
                analyzer.ctx.has_dapi = ch_cfg.get("has_dapi", False)
                analyzer.ctx.dapi_channel_idx = ch_cfg.get("dapi_idx")
                analyzer.ctx.has_granules = ch_cfg.get("has_granules", False)
                analyzer.ctx.granule_channel_idx = ch_cfg.get("granule_idx")

                analyzer.ctx.selected_channels = [analyzer.ctx.cell_channel_idx]
                if analyzer.ctx.has_dapi and analyzer.ctx.dapi_channel_idx is not None:
                    analyzer.ctx.selected_channels.append(analyzer.ctx.dapi_channel_idx)
                if analyzer.ctx.has_granules and analyzer.ctx.granule_channel_idx is not None:
                    analyzer.ctx.selected_channels.append(analyzer.ctx.granule_channel_idx)

                analyzer.ctx.selected_channels = sorted(list(set(analyzer.ctx.selected_channels)))
                analyzer.ctx._channels_initialized = True
            else:
                analyzer.ctx._channels_initialized = False

            if not analyzer.ctx._channels_initialized:
                self.root.after(0, lambda: self.append_log(self._t("log.asking_channel_choice")))
                analyzer.optics.determine_channels(files[0])
                analyzer.ctx._channels_initialized = True

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
