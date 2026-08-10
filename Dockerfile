FROM python:3.10-slim

ENV OPENSSL_FORCE_FIPS_MODE=0 \
    OPENSSL_FIPS=0

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN echo '#define _GNU_SOURCE\n\
#include <stdio.h>\n\
#include <string.h>\n\
#include <dlfcn.h>\n\
#include <stdarg.h>\n\
typedef FILE *(*fp_t)(const char *, const char *);\n\
typedef int (*op_t)(const char *, int, ...);\n\
FILE *fopen(const char *p, const char *m) {\n\
    fp_t o = (fp_t)dlsym(RTLD_NEXT, "fopen");\n\
    return (p && strcmp(p, "/proc/sys/crypto/fips_enabled") == 0) ? o("/dev/null", m) : o(p, m);\n\
}\n\
FILE *fopen64(const char *p, const char *m) {\n\
    fp_t o = (fp_t)dlsym(RTLD_NEXT, "fopen64");\n\
    if (!o) o = (fp_t)dlsym(RTLD_NEXT, "fopen");\n\
    return (p && strcmp(p, "/proc/sys/crypto/fips_enabled") == 0) ? o("/dev/null", m) : o(p, m);\n\
}\n\
int open(const char *p, int f, ...) {\n\
    op_t o = (op_t)dlsym(RTLD_NEXT, "open");\n\
    if (p && strcmp(p, "/proc/sys/crypto/fips_enabled") == 0) return o("/dev/null", f);\n\
    va_list a; va_start(a, f); int m = va_arg(a, int); va_end(a);\n\
    return o(p, f, m);\n\
}\n\
int open64(const char *p, int f, ...) {\n\
    op_t o = (op_t)dlsym(RTLD_NEXT, "open64");\n\
    if (!o) o = (op_t)dlsym(RTLD_NEXT, "open");\n\
    if (p && strcmp(p, "/proc/sys/crypto/fips_enabled") == 0) return o("/dev/null", f);\n\
    va_list a; va_start(a, f); int m = va_arg(a, int); va_end(a);\n\
    return o(p, f, m);\n\
}' > /tmp/nofips.c && \
    gcc -shared -fPIC /tmp/nofips.c -o /usr/local/lib/libnofips.so -ldl && \
    rm /tmp/nofips.c

ENV LD_PRELOAD=/usr/local/lib/libnofips.so

WORKDIR /app

RUN pip install --no-cache-dir setuptools wheel Cython pysam numpy pandas scipy

COPY . /app

RUN pip install --no-cache-dir --no-build-isolation .

ENTRYPOINT ["indelinside"]
