/* fnv1a64: независимый пересчёт FNV-1a 64 для проверки CBT1-файлов.
 * Использование: fnv1a64 FILE OFFSET LENGTH  ->  hex-дайджест байт
 * [OFFSET, OFFSET+LENGTH). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s FILE OFFSET LENGTH\n", argv[0]);
        return 3;
    }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror("fopen"); return 3; }
    long long off = atoll(argv[2]);
    long long len = atoll(argv[3]);
    if (fseek(f, off, SEEK_SET) != 0) { perror("fseek"); return 3; }

    uint64_t h = 0xCBF29CE484222325ULL;
    char buf[1 << 20];
    long long left = len;
    while (left > 0) {
        size_t want = left < (long long)sizeof(buf)
                          ? (size_t)left : sizeof(buf);
        size_t got = fread(buf, 1, want, f);
        if (got == 0) { fprintf(stderr, "short read\n"); return 3; }
        for (size_t i = 0; i < got; ++i) {
            h ^= (unsigned char)buf[i];
            h *= 0x100000001B3ULL;
        }
        left -= (long long)got;
    }
    fclose(f);
    printf("%016llx\n", (unsigned long long)h);
    return 0;
}
