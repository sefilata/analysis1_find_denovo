CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall

BIN_DIR := exec

TARGETS := \
	$(BIN_DIR)/03_validate_count.out \
	$(BIN_DIR)/03-2_validate_count_GC.out

all: $(TARGETS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/03_validate_count.out: 03_validate_count.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BIN_DIR)/03-2_validate_count_GC.out: 03-2_count_GC.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS)

.PHONY: all clean
