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

### Canonical K3 storage and candidate-pair collapse

`--glv-orbit --canonical-k3` stores one canonical row per C6 orbit instead of
one row per table scalar. The stored scalar is the orbit minimum

```text
c = min over r in {0,1,2} of min(lambda^r * k mod n, n - lambda^r * k mod n),
```

which is constant on the orbit and distinct across orbits. `K3 = x^3` and the
`(d01, d12)` signature depend on x-coordinates only, so the canonical
representative keeps the same bucket key, and its signature is read off the
phi-rotated x-coordinate `beta^r * x` of the row point (negation does not
change x), so canonicalization costs no extra scalar multiplication.
With `--table-scalars` the build is a direct canonical pass: the input list
is canonicalized and deduped BEFORE any point is computed (two order
multiplications per scalar, no EC work), so exactly one point per orbit is
multiplied and the build cost and peak memory follow the stored table, not
the raw input list. `merged rows` counts input rows dropped
before the first scalar multiplication and `table rows stored` shows the
final size. On the continuous window table `[0, m)` the merge is a no-op
(partners `lambda*j ~ 2^128` never co-occur there); it becomes active with
`--table-scalars`.

Every `ComputePublicKey` (build rows, anchor/step/shift points, exact
candidate checks) and every real point addition is counted; the output
reports `EC scalar multiplications`, `EC point additions` and
`total EC operations`, so builds and search modes compare on the same unit.

The phase-signature path now performs candidate-pair level dedup inside a
bucket: every row routes to its pair `{+o, -o}`, identical pairs are skipped
(`signature pair collapses`), and differing pairs are still verified
defensively (`signature pair mismatches`, expected zero — all rows of a K3
bucket lie in one C6 orbit and route to the same pair).

`--table-scalars FILE` builds the baby table from an explicit list of hex
scalars (one per line, `#` comments allowed, zero rejected) instead of the
continuous window `[0, m)`. The giant plan strides by the stored row count,
so targets `start + i*m + s` with `s` in the orbit set of the table remain
recoverable. This is the research mode in which K3 buckets genuinely contain
orbit partners — the full-range representative table anticipated by the
routing-economy measurements — and it instantiates the storage/recovery
composition: canonical storage gives one row per orbit, partner-expansion
(`--partner-table`) restores the target from the routed pair.

```bash
# table with a real 3-row bucket (j3, lambda*j3, lambda^2*j3 plus background)
./eng_bsgs --target PUBKEY --start 0x4000000000000000000000000000000000 \
  --size 1048576 --segment 1048576 --glv-orbit --phase-signature \
  --canonical-k3 --partner-table --table-scalars rows.txt --threads 4
```

Verified by `test_canonical_k3.py` (negative-control collapse counting,
partner recovery, canonical merge 3 rows -> 1, in-window and continuous-window
regressions) and `economy_engine_sweep.py` (24 targets: baseline six-branch
3.04 exact EC checks per hit vs 1.33 routed, ~2.28x; canonical storage halves
the table with the same pair economy).

### Representatives file for long runs

`make_representatives_file.py` generates persistent canonical
representative files (`representatives_canonical_H*.txt`) for
`--table-scalars`: one canonical orbit representative per row, heads
`j = 1..H`, deterministic and prefix-nested across H. For heads `j <= H`
the canonical representative equals `j` itself (the lambda-rotated
partners are of order `n`), so the file's core is the contiguous tiling
`i*stride + [1, H]`; the full range comes from orbit partners
`+-lambda^a * j` at every giant position via `--partner-table`:

```bash
./eng_bsgs --target PUB --start ... --size CHUNK --segment CHUNK \
  --glv-orbit --phase-signature --partner-table \
  --table-scalars ../representatives_canonical_H400000.txt
```

Two operational facts validated by the long cover run
(`long_cover_validation.py`, PROGRESS_MEMO.md):

- **Set `--segment` explicitly.** Without `--baby-rows` the default
  `segment_size` is 100000; a `--size 5e11` range is then split into
  5M segments, each paying its own anchor/step scalar multiplications
  (the giant walk itself does none) - measured 17.8k giant steps/second
  instead of ~330k, an ~18x slowdown. With `--segment CHUNK` (one
  segment per cover chunk) the full ~330k steps/second is restored.
- **Stride under `--canonical-k3`.** The canonicalization pre-pass runs
  before `m = table_scalars.size()`, so the stride is the POST-merge row
  count (1.2M input rows merge to 400k, `merged rows: 800000`),
  not the file's line count. The merged mode therefore has the same
  stride, giant-step count and full coverage as the canonical file;
  the plain 3H-row file without `--canonical-k3` strides by 3H and
  covers only `[1, H]` after each giant position (holes in `(H, 3H)`),
  which the negative controls of the long cover run demonstrate.

The giant stride equals the stored row count. The H-curve measurements
(50k/100k/200k/400k), the long cover validation (13/13 scenarios, 27
SIGKILL/resume cycles, hole-free and double-count-free chunk accounting)
and the black-box round-trip verification are covered by
`test_representatives_file.py`, `long_cover_validation.py` and
PROGRESS_MEMO.md; file checksums are pinned in
`../representatives_files.json`.

Format roles of the representatives artifacts. Today the engine reads the
`.txt` file directly through its own parser (`--table-scalars`: one hex
scalar per line, `#` comments, `0x`/hex or decimal accepted, zero
rejected). There is no separate binary engine format yet, so a single
`.txt` artifact currently plays two roles:

1. diagnostic export — human-readable, diffable, checksummed in
   `representatives_files.json`;
2. engine input — parsed per line at startup.

For H <= 400k this is fine (400k rows ≈ 3.1 MB of text, parsing is a
sub-second fraction of the build). For H = 4M the `std::vector<Int>`
input list alone costs about 160 MB plus parse time, and merge-mode
expansion (`--canonical-k3` on an expanded list) triples that, so the
agreed upgrade path for H=4M is a dedicated binary engine format
(`canonical_version` header, sorted fixed-width rows, checksum,
mmap-friendly read) with the `.txt` kept as the diagnostic export. Until
that format exists, H = 4M runs must not rely on the raw `.txt` input
path.

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

### Checkpoint/resume для длинных --cover прогонов

Длинные прогоны проектируются под убивающую среду с первого дня. За
`--cover-start/--cover-size` закреплены три опции:

```bash
./eng_bsgs --target PUBKEY --start 0x4000000000000000000000000000000000 \
  --size 1000000000000 --baby-rows 1000000 \
  --cover-start 0x4000000000000000000000000000000000 \
  --cover-size 3000000000000 \
  --glv-orbit --phase-signature --partner-table \
  --table-scalars ../representatives_canonical_H400000.txt \
  --checkpoint-file /scratch/cover.cp --checkpoint-steps 1 \
  --threads 16
```

- `--checkpoint-file PATH` — атомарный checkpoint прогресса: каждые
  `--checkpoint-steps N` завершённых внешних chunk (N>=1, по умолчанию 1 —
  после каждого chunk, т.е. хотя бы раз за ~12 с при H=400k) запись в
  `PATH.tmp` → `fsync` → `rename`. Гранулярность — внешний chunk: это один
  гигант-проход шириной `--size`; прерванный chunk перепрокрывается
  целиком, поэтому двойных дыр и двойного счёта нет.
- В файле хранятся: plan fingerprint (FNV-1a от target/start/size/cover/
  режимов/файла таблицы), позиция следующего chunk, число завершённых
  chunk, накопительные EC/шаговые счётчики, флаг `found` и (при успехе)
  восстановленный ключ.
- `--resume` продолжает план с позиции checkpoint. Несовпадение
  fingerprint (другой target/размер/cover/режимов) отклоняется с кодом 2.
  Поисковые счётчики восстанавливаются, поэтому финальный вывод
  кумулятивен между сессиями (build-строка — только текущая сессия: таблица
  строится заново в каждой).
- Если checkpoint сообщает `found: 1`, `--resume` печатает найденный ключ
  из checkpoint и не ищет заново; правило «прогон не считается успешным
  без ≥1 найденной цели» контролирует обходной harness (контрольные цели
  засеиваются в покрытый диапазон, см. `test_cover_checkpoint.py`).
- Fallback-проход (zero shift) после исчерпания cover тоже попадает в
  checkpoint.

Верификация: `test_cover_checkpoint.py` (black-box) — SIGKILL после
checkpoint chunk 0 → `--resume` находит контрольную цель без переделки
chunk 0 (`automatic chunks attempted: 1`) → повторный `--resume` печатает
результат без поиска → `--resume` чужого плана отклоняется по
fingerprint. ALL PASS.

### Split EC-accounting и TTS в выводе

В выводе счётчики EC разделены по фазам:

```text
EC scalar multiplications: <total, кумулятивен при --resume>
EC scalar multiplications (build): <fresh build этой сессии>
EC scalar multiplications (giant search, cumulative): <total - build - hit>
EC scalar multiplications (exact candidate checks): <hit-путь: ec_mult>
EC point additions / (build) / (giant search, cumulative)
total EC operations: <total>
```

При успехе дополнительно печатаются:

```text
time to solution seconds: <с момента старта, включает build>
exact candidate checks per hit: <ratio exact/search hits> (exact E / search hits H)
```

### Счётчики по фазам: merged rows / insert collisions / search hits
### (30.10, семантика зафиксирована)

Раньше один и тот же счётчик за сессию менял имя и смысл (build hits ->
canonical rows merged -> «это вставки»), и в merge-логе
`build hits: 800000` читалось как «движок нашёл 800k целей до начала
поиска». Теперь имена говорят, что именно считается:

```text
merged rows:       <входные строки, схлопнутые в существующую строку
                   при СБОРКЕ таблицы (орбита-коллапс: равный канонический
                   скаляр или равный K3-ключ). Ожидаемо по построению,
                   НЕ поисковое событие: 2 на голову на 3-строчном
                   expanded-файле, 0 на canonical-файле и непрерывном
                   окне. Session-pure, --resume не восстанавливает>
insert collisions: <строки с равным K3-ключом, но другим сохранённым
                   скаляром - настоящая межорбитальная коллизия. На полных
                   C6-орбитах обязана быть 0 (равный K3-ключ = ровно одна
                   орбита); ненулевое значение - сломанный инвариант:
                   строка отбрасывается, stderr предупреждает. Детектируется
                   только в --canonical-k3, где сохранённый скаляр -
                   минимум орбиты>
search hits:       <бакет-совпадения гигант-скана, один на совпавшую
                   гигант-позицию независимо от числа строк в бакете;
                   восстанавливается --resume, кумулятивен>
invariant hits:    <legacy row-level счётчик: строки совпавших бакетов,
                   включая пропущенные дедупом>
```

`merged rows == 0` на canonical-файле - инвариант каноничности (файл без
дублей орбит); `merged rows == 2 * heads` в merge-режиме (1.2M входных
строк -> 400k орбит). `exact candidate checks per hit` - ratio
exact/search hits, переименованием не задет. `duplicate baby keys`
остаётся метрикой качества таблицы (raw vs unique), а не счётчиком
событий. `giant steps/second` - сессионно-честный: из кумулятивного
счётчика вычитаются шаги, восстановленные из checkpoint.

### Прогресс-тикер и assert стрида (30.10)

`--progress-interval N` (секунды, default 60, 0 = off): поток движка
раз в N секунд поиска печатает

```text
progress: elapsed 60.0s | giant steps 19628898 | steps/s 326933.3 | rss kB 81760
```

шаги/с - сессионные (минус restored offset), RSS - текущий VmRSS. Строки
предназначены для внешних логов стабильности длинных прогонов.

Assert стрида: после сборки из файла представителей движок проверяет
`stored + canonical_rows_merged == raw_file_rows` (иначе exit 2) и печатает
`giant stride (post-merge rows): <N>` - стрид гигант-плана равен строкам
ПОСЛЕ merge (для canonical-файла H400000 это 400000). Нюанс чтения лога:
`giant steps` кумулятивен по воркерам; при находке на границе сплита
воркеров счётчик показывает шаги другого воркера (единицы-десятки), а не
позицию находки - это не аномалия.

Каждая точная проверка кандидата (`ec_mult(candidate)`, `ec_mult(recovered)`
в шести-ветвевом, маршрутизованном, partner и plain путях) учитывается в
hit-счётчике; гигант-строка = search − hit. `cache-misses` в движке не
собирается: perf в среде недоступен (см. PROGRESS_MEMO.md «perf: не
измерено»); если perf станет доступен, счётчик снимается обходным
процессом на throughput-фазе.

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
