CC=gcc
LD=gcc

SRC_PATH=./src/
OUTPUT=sar_image_formation_tool

TOOL_BIN_PATH=./bin/
TOOL_OBJECTS=tool/main.o

OBJECTS=$(SRC_PATH)/sar_image_formation_tool.o
HEADERS=

CINCLUDE_PATH= -I$(SRC_PATH)
CFLAGS=-Wall -fPIC $(CINCLUDE_PATH)
LDFLAGS=-Wl,--no-as-needed -lm -lers_raw_parser -lsarif -lmatio

.PHONY: all release clean

all: CFLAGS += -g
all: $(OUTPUT)

release: CFLAGS += -O2
release: $(OUTPUT)

clean:
	- $(RM) $(OUTPUT) $(OBJECTS) *~ $(SRC_PATH)/*~

$(OUTPUT): $(OUTPUT_PATH) $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.c $(HEADERS)
	$(CC) -c $(CFLAGS) -o "$@" "$<"
