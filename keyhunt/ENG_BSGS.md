# Experimental engineering BSGS

`eng_bsgs` is a separate benchmark built on the existing `keyhunt/secp256k1`
implementation. It deliberately does not modify the production `keyhunt`
search modes or their on-disk Bloom/table format.

The program covers a finite interval completely. It partitions the interval
into segments, builds a local baby table for one segment, searches the segment
with giant steps, and then reuses the table for the next segment. This exposes
the memory/time trade-off measured in the Python experiments.

Build it with:

```bash
cd keyhunt
make eng_bsgs
```

Run a bounded test with a known compressed public key:

```bash
./eng_bsgs \
  --target 02...64_hex_x... \
  --start 1 \
  --size 300000 \
  --segment 10000 \
  --threads 4
```

The target format is SEC1 hex: compressed keys begin with `02` or `03`, and
uncompressed keys begin with `04`. `--start` accepts decimal or `0x...`
hexadecimal text. `--size` is the exact finite interval width; the searched
interval is `[start, start + size)`.

Output includes:

- number of segments and the shared read-only baby table;
- one-time baby-table construction size;
- giant-step count and giant-steps/second;
- recovered scalar, segment and local offset when a match is found.

The `recovered private key` line prints the same scalar zero-padded to 64 hex
characters, so it can be copied directly as a 32-byte secp256k1 private key.

The benchmark stops early on the first exact point match. A false result means
that the supplied public key was not found in the specified finite interval.
It does not claim anything about the rest of the secp256k1 key space.

### Lambda-sector priority

The experimental Lambda router can order the numeric slices of each local
window using the three exact `w1.py` Lambda values:

```bash
./eng_bsgs --target PUBKEY --start 1 --size 1000000 --segment 10000 \
  --lambda-priority --lambda-sectors 16
```

The three sectors are searched first. Every other sector is then searched as a
full fallback, so this option changes only traversal order and cannot remove a
candidate. The output reports Lambda calculations and the number of priority
and fallback slices separately; those counters are experimental statistics,
not a mathematical speed guarantee. `--glv-orbit` can be used together with
this mode.

### Phase-signature GLV branch routing

`--glv-orbit --phase-signature` refines the K3=x^3 quotient with the ordered
Lambda signature `(d01, d12)`, where `d01 = Lambda(phi(P)) - Lambda(P)` and
`d12 = Lambda(phi^2(P)) - Lambda(phi(P))`. The signature is invariant under
`P -> -P` and cyclically permuted by `phi`, so a K3 bucket match can be routed
to the single branch pair `{+lambda^r * j, -lambda^r * j}` instead of all six
exact branches. Degenerate rows fall back to the full six-branch check, so
coverage never depends on the router.

```bash
./eng_bsgs --target PUBKEY --start 0x4000000000000000000000000000000000 \
  --size 100000 --segment 100000 --glv-orbit --phase-signature --threads 4
```

The output adds `signature routes` and `signature anomalies (6-branch
fallback)` counters, plus the real `unique baby rows built` /
`duplicate baby keys` table statistics. Measured on a bounded 135-bit window,
routed recovery is correct and route anomalies are zero; the measured exact
check count does not improve there because the `r=0` branch is tried first in
both modes. The 6-to-2 candidate reduction is verified structurally by
`phase_signature_135bit_experiment.py` (orbit mode: 1500/1500 correct routes,
1500/1500 candidate pairs contain the target offset). The table build costs
about three Lambda computations per row in this mode.

### Partner-table recovery of out-of-window orbit branches

`--glv-orbit --partner-table` widens the recovery contract: orbit-branch
candidates that fall outside the searched interval are no longer skipped by
the interval pre-filter but are exact-verified against the target with a
full EC check. This is safe because a K3 hit never fires spuriously — the
roots of `X^3 - c` over F_p are exactly `x, beta*x, beta^2*x`, so a bucket
match means the giant point IS `±lambda^a * jG` for a table row `j`, and a
candidate `i*m ± lambda^a * j` that passes `ComputePublicKey` equals the
target scalar. The mode recovers scalars of the form

```text
k = range_start + i*m ± lambda^r * j (mod n),  j < baby rows,
```

i.e. orbit partners of window scalars, which the default mode silently
drops. The cover-range check is bypassed in this mode (the EC verification
against the original target is the ultimate guard); partner-type recoveries
are counted in the new `partner recoveries` output line. Verified by
`test_partner_table.py` (black-box: negative control without the flag,
six-branch recovery, signature-routed recovery, in-window regression).

```bash
./eng_bsgs --target PUBKEY --start 0x4000000000000000000000000000000000 \
  --size 1048576 --segment 1048576 --glv-orbit --phase-signature \
  --partner-table --threads 4
```

### Automatic external coverage

`--cover-start` and `--cover-size` describe an external scalar interval using
the full `Int` representation (so values such as `2^134` are accepted). The
local BSGS window is translated over that interval one chunk at a time:

```bash
./eng_bsgs --target PUBKEY --start 1 --size 100000 \
  --segment 10000 --cover-start 0x4000000000000000000000000000000000 \
  --cover-size 0x4000000000000000000000000000000000 \
  --lambda-priority --lambda-sectors 16
```

For a source chunk beginning at `s`, the program uses
`t = start - s (mod n)` and searches `Q+tG`. A match is checked against both
the shifted target and the original target, then returned as `k = k_t - t
(mod n)`. The final short chunk uses its actual width. The original target is
always searched last as a fallback. Manual `--shift` arguments and automatic
coverage are rejected together to keep the search schedule unambiguous.

The focused regression harness checks the exact Lambda values, priority-only
search, the last short automatic chunk, GLV combination, and argument validation:

```bash
python3 test_lambda_priority.py
```

Для полного сценария «сгенерировать пару, покрыть внешний сектор и сдвинуть
цель»:

```bash
python3 test_lambda_cover.py
```

Этот тест генерирует случайный scalar в 135-битном диапазоне, помещает его в
последний неполный чанк внешнего интервала и передаёт программе только
публичный ключ. Включены одновременно Lambda-приоритет и `--glv-orbit`.

Чтобы явно разделить «размер таблицы» и «покрываемый диапазон», можно задать
число baby-строк `B`. Тогда один сегмент имеет ширину примерно `B²`:

```bash
./eng_bsgs \
  --target 03...64_hex_x... \
  --start 0x4000000000000000000000000000000000 \
  --size 3000000000000 \
  --baby-rows 1000000 --threads 4
```

Здесь таблица содержит около миллиона baby-строк на worker, а покрываемое окно
составляет три триллиона кандидатов. Это уже настоящий BSGS trade-off: работа
порядка `sqrt(3e12)` giant-шагов в одном сегменте. `--segment N` остаётся
альтернативой, если нужно задать ширину сегмента напрямую.

Если ограничить resident-таблицу примерно 300 000 baby-строк:

```bash
./eng_bsgs \
  --target YOUR_COMPRESSED_PUBLIC_KEY \
  --start 0x4000000000000000000000000000000000 \
  --size 3000000000000 --baby-rows 300000 --threads 4
```

Один сегмент тогда имеет ширину примерно `300000² = 9e10`, поэтому всё окно
из `3e12` значений покрывается примерно 34 сегментами. В измерении с целью,
находящейся около `2.5e12` от начала окна, программа действительно вывела
`segments: 34`, `baby rows (shared): 300000`, `recovery correct: true`.
Это и есть требуемая схема «таблица 300 тысяч, покрытие триллионов»: память
ограничена одной таблицей, а цена — повторный BSGS-проход по сегментам.

Для проверки именно нашей GLV-гипотезы добавлен флаг `--glv-orbit`. Он заменяет
точный ключ таблицы на `K3=x^3`, а после совпадения проверяет шесть точных
скалярных ветвей `±j`, `±lambda*j`, `±lambda²*j`. Это не «магическое» удаление
кандидатов: таблица хранит только orbit quotient, а точность подтверждается
вычислением public point для каждого кандидата.

Black-box запуск с GLV:

```bash
python3 test_eng_bsgs.py --seed 20260925 --glv-orbit
```

### Что даёт «размножение» публичной точки

Для цели `Q = kG` можно вычислить `tQ = (tk)G` и искать скаляр `s = tk`,
после чего вернуть `k = t⁻¹s mod n`. Это полезно как способ построить
несколько связанных потоков поиска, но не добавляет независимых бит информации:
при известном `t` это просто обратимая перенумерация скаляров. В частности,
для короткого окна `[a, a+w)` новое окно имеет вид `[ta, ta+tw)` (с учётом
модуля `n`) и обычно не становится меньше.

Режим `--glv-orbit` проверяет более содержательный вариант: для каждого
giant-point используется инвариант `K3=x^3`, общий для шести точек
`±Q, ±φ(Q), ±φ²(Q)`, а затем выполняется точная проверка допустимых ветвей.
Это может уменьшить постоянный размер таблицы на полном диапазоне, но не
делает giant-проход субкорневым. На коротком bounded-окне ветви `λk` и
`λ²k` обычно лежат вне окна, поэтому выигрыш часто отсутствует или становится
отрицательным из-за стоимости exact-проверок.

Важная деталь реализации: одинаковый `K3` не удаляется из baby-таблицы.
Все строки одного bucket проверяются точно; это исключает ложное ускорение за
счёт пропуска кандидатов.

Чтобы исключить подгонку public key под заранее известный scalar, используйте
black-box harness. Он генерирует новый случайный scalar внутри окна, передаёт
в `eng_bsgs` только public key/start/size, а затем сравнивает напечатанный
64-значный результат с контрольным scalar внутри самого тестового процесса:

```bash
python3 test_eng_bsgs.py --seed 20260925
```

Сам `eng_bsgs` не получает private key и не добавляет target в baby-table:
таблица строится только из `j*G`, а совпадение проверяется отдельным
вычислением `candidate*G == target`.

The baby table is built once and shared read-only by all worker threads. Earlier
versions built an identical private copy inside every worker; that multiplied
RAM usage by `--threads` without increasing coverage. Sharing the table changes
memory use, not the BSGS work or the mathematical coverage.

This is an engineering benchmark, not a new sub-root discrete-log algorithm.
Reducing the resident table only moves work into rebuilding/searching more
segments. The original `keyhunt` binary and its normal `bsgs` mode remain
unchanged.

## Команда настоящего keyhunt для полного 135-битного окна

Для полного окна `[2^134, 2^135)` используй штатный бинарник `keyhunt`, а не
`eng_bsgs`. Создай файл с одним своим public key:

```bash
printf '%s\n' 'YOUR_COMPRESSED_PUBLIC_KEY' > target.pub
```

Затем запусти последовательный BSGS:

```bash
./keyhunt \
  -m bsgs -f target.pub \
  -b 135 -B sequential \
  -n 0x100000000000 -k 512 \
  -t 32 -q -s 10 -S
```

Здесь `-b 135` задаёт именно диапазон `[2^134,2^135)`, `-n` — ширину одного
BSGS-блока, а `-k 512` увеличивает resident-таблицу по схеме из README. `-S`
сохраняет Bloom/table-файлы для повторного запуска; перед использованием нужно
проверить свободные RAM и диск. `-R` здесь не нужен: он переключает поиск в
случайный режим и не даёт последовательного покрытия окна.

Строки `Pkeys/s` в README — это внутренняя BSGS-метрика эквивалентного
покрытия, а не буквальный перебор каждого scalar. Реальное время для полного
окна определяется статистикой запущенного процесса и остаётся огромным даже
при высокой метрике.

## Приоритет сдвинутых целей

`eng_bsgs` принимает повторяемый параметр `--shift`:

```bash
./eng_bsgs --target PUBKEY --start 0x1000 --size 0x100 \
  --shift=-0x4000 --shift=-0x3000
```

Для каждого `t` выполняется поиск цели `Q+tG`. Если найден локальный скаляр
`k_t`, программа восстанавливает и проверяет исходный:

```text
k = k_t - t (mod n)
```

Исходная цель всегда проверяется последней как fallback, поэтому добавление
сдвигов не удаляет существующее покрытие. Каждый сдвиг является отдельным
BSGS-проходом и может увеличить работу. Режим предназначен для эксперимента
с приоритетом и покрытием, а не доказывает, что фиксированного числа сдвигов
достаточно для гарантированного восстановления.

Для воспроизводимой проверки 135-битной арифметики есть harness:

```bash
python3 test_135bit_shift.py --seed 20260926
```

Он генерирует случайный `k` в `[2^134, 2^135)`, строит compressed public key,
переводит его известным тестовым сдвигом в маленькое окно и проверяет возврат
исходного `k`. Знание этого сдвига относится только к тестовому harness; это не
слепое восстановление произвольного 135-битного ключа.
