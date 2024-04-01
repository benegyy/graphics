rasterizer:
# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -std=c11

# Source files

SRCS = Vec4.cpp Vec3.cpp my_line.cpp Triangle.cpp Translation.cpp Color.cpp Matrix4.cpp Helpers.cpp my_transformation.cpp Mesh.cpp Rotation.cpp Scaling.cpp Camera.cpp tinyxml2.cpp Scene.cpp Main.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
TARGET = rasterizer

# Default rule
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -lstdc++ -lm

# Rule to compile source files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
