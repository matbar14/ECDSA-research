#!/usr/bin/env python3
"""Числа для ответа на вопрос: можно ли хранить BSGS-таблицы в кубитах.

1) Holevo: n кубитов несут максимум n бит ДОСТУПНОЙ классической информации
   (фон Нейман-энтропия n кубитов <= n). Значит таблица из B байт требует
   >= 8*B кубитов в принципе, без всякой структуры запросов.
2) QRAM (бинарное дерево GLM'08): M записей по L бит -> O(M*L) кубитов
   (неструктурированные данные не сжимаются), запрос в суперпозиции за O(log M).
3) Сравнение с нашими таблицами из h_curve_summary.json и с огибающей
   statevector-симуляции (32 GiB -> 31 кубит).
4) secp256k1: классическая таблица против квантовой схемы Шора
   (Roetteler et al. 2017: ~2330 логических кубитов на 256-битный ECDLP).
"""
import json

# --- наши таблицы из H-кривой ---
hs = json.load(open("/home/z/my-project/ECDSA/results/h_curve_summary.json"))
print("== Наши BSGS-таблицы (H-кривая) ==")
for pt in hs["points"]:
    h = pt["heads"]
    c = pt["throughput"]["canonical H rows"]
    rows, bpr, rss = c["table rows stored"], c["table bytes per row"], c["peak rss mb"]
    bits = rows * bpr * 8
    print(f"H={h:>7}: {rows:>7} строк x {bpr:5.1f} B = {rows*bpr/2**20:6.1f} MB "
          f"= {bits/2**20:7.1f} Mбит -> минимум {bits/8/2**20:6.1f} млн кубитов "
          f"(Holevo), RSS {rss} MB")

# --- огибающая 32 GB ---
env = json.load(open("/home/z/my-project/ECDSA/results/qip_statevector_envelope.json"))
qt = env["qubits max target (GiB=32)"]
print(f"\n== Огибающая statevector: 32 GiB -> {qt} кубитов "
      f"({16*2**qt/2**30:.0f} GiB амплитуд) ==")
print(f"   31 кубит несёт максимум {qt} бит доступной информации (Holevo)")

# --- таблица H=400k против регистра симуляции ---
canon = hs["points"][3]["throughput"]["canonical H rows"]
table_bits = canon["table rows stored"] * canon["table bytes per row"] * 8
sv_bits_capacity = qt  # Holevo
print(f"\nТаблица canonical H=400k: {table_bits/2**20:.0f} Мбит")
print(f"Регистр нашей симуляции (31 кубит): {sv_bits_capacity} бит доступно")
print(f"Дефицит: x{table_bits/sv_bits_capacity/2**20:.1f} млн")

# --- secp256k1 в масштабе ---
print("\n== secp256k1 (n ~ 2^256) ==")
for name, rows in (("наш режим H=400k", 400_000),
                   ("H=10^9 (полноразмерный BSGS)", 10**9),
                   ("baby steps 2^128 (чистый BSGS)", 2**128)):
    b = rows * 340  # байт на строку по нашей кривой
    print(f"{name}: таблица {b/2**30:.1e} GB, минимум {b*8/2**40:.2e} млрд кубитов "
          f"(Holevo)")
print("Схема Шора на 256-битный ECDLP (Roetteler et al. 2017): "
      "~2330 логических кубитов, таблица не нужна вовсе")
print(f"Разрыв с железом: 10^3 физических кубитов сегодня против "
      f"10^9+ логических для таблицы, но 2.3e3 для Шора")
