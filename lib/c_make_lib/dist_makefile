CC=g++
CFLAGS=-O3 -Wno-write-strings
SOURCES := $(shell find . -type f -name *.c)
OBJECTS := $(SOURCES:.c=.o)
DEPS := $(OBJECTS:.o=.deps)
 
t_coffee: $(OBJECTS)
	@echo " Linking..."; $(CC) $^ -o t_coffee -lm
 
%.o: %.c
	@echo " CC $<"; $(CC) $(CFLAGS) -I. -MD -MF $(@:.o=.deps) -c -o $@ $<
 
-include $(DEPS)
 
all: t_coffee TMalign
 
TMalign:
	$(FCC) TMalign.f -o TMalign 
 
.PHONY: clean

clean:
	@echo " Cleaning..." 
	`find . -name *.o | xargs rm -rf {}`