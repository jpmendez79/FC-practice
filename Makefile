CXX = g++
CXXFLAGS = -g -Wall -lMinuit `root-config --cflags --libs `

# List of source files (add more if needed)
SOURCES = matrix.cpp

# Name of the executable
TARGET = matrix

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

# Clean rule
clean:
	rm -f $(TARGET)
