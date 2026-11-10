#!/bin/bash
# H-кривая масштабирования: 50k / 200k / 400k голов на честной базе
# (прямой canonical-проход + total EC-ops), 2 потока = условия базы H=100k.
# Использование: run_h_curve.sh <heads>
set -e
HEADS="$1"
cd /home/z/my-project/ECDSA/keyhunt
python3 table_representative_sweep.py \
    --heads "$HEADS" \
    --threads 2 \
    --output "../h_curve_H${HEADS}.json"
echo "[done] H=${HEADS}"
