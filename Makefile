CC=clang++
CCFLAGS=-std=c++17 -Wall -Wextra -isystem stb -I . # flags for g++: -lm -lstdc++fs

MAIN=transcoder
DEPS=*.h stb/*.h

.PHONY: clean force debug release bear

$(MAIN): $(MAIN).o
	$(CC) $(CCFLAGS) $^ -o $@

$(MAIN).o: $(MAIN).cpp $(DEPS)
	$(CC) -c $(CCFLAGS) $< -o $@

debug: CCFLAGS += -g -O0
debug: clean $(MAIN)

release: CCFLAGS += -O3 -DNDEBUG
release: clean $(MAIN)

quantizer: quantizer.cpp
	$(CC) $(CCFLAGS) -O3 quantizer.cpp -o quantizer

clean:
	rm -f $(MAIN).o $(MAIN)

bear: clean
	bear make

force: clean $(MAIN)
