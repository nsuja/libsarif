CC=gcc
LD=gcc
INSTALL=install

INSTALL_HEADERS_PATH=/usr/include/
INSTALL_LIB_PATH=/usr/lib/

SRC_PATH=./src/
OUTPUT_PATH=./lib/
OUTPUT=libsarif.so

OBJECTS=src/sarif.o
HEADERS=src/sarif.h

CINCLUDE_PATH= -I$(SRC_PATH)
CFLAGS=-g -Wall -fPIC $(CINCLUDE_PATH) `pkg-config --cflags fftw3`
LDFLAGS=-shared -Wl,--no-as-needed,-soname,$(OUTPUT) -lm  `pkg-config --libs fftw3`

all: $(OUTPUT)

.PHONY: tool clean install

tool:
	make -C tool/ all

clean:
	- $(RM) $(OUTPUT) $(OBJECTS) *~ $(SRC_PATH)/*~
	- make -C tool/ clean

install:
	- $(INSTALL) "$(OUTPUT_PATH)/$(OUTPUT)" $(INSTALL_LIB_PATH)
	- $(INSTALL) $(HEADERS) $(INSTALL_HEADERS_PATH)

$(OUTPUT): $(OUTPUT_PATH) $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o "$(OUTPUT_PATH)/$@"

%.o: %.c $(HEADERS)
	$(CC) -c $(CFLAGS) -o "$@" "$<"

$(OUTPUT_PATH):
	- mkdir $(OUTPUT_PATH)
