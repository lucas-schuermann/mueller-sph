#######################################################################################################

# Mac OS X
INCLUDE_PATH      = -I/usr/local/include/
LIBRARY_PATH      = -L/usr/local/lib/
OPENGL_LIBS       = -framework OpenGL -framework GLUT

# # Linux
#INCLUDE_PATH      =
#LIBRARY_PATH      =
#OPENGL_LIBS       = -lglut -lGL -lX11

# # Windows / Cygwin
# INCLUDE_PATH      = -I/usr/include/opengl
# LIBRARY_PATH      = -L/usr/lib/w32api
# OPENGL_LIBS       = -lglut32 -lopengl32

#######################################################################################################

TARGET = sph
CC = g++
LD = g++
CFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated -pedantic -Wno-vla-extension $(INCLUDE_PATH) -I./include -I./src -DNDEBUG
LFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated -Werror -pedantic $(LIBRARY_PATH) -DNDEBUG
LIBS = $(OPENGL_LIBS)

OBJS = obj/main.o

default: $(TARGET)

all: clean $(TARGET)

$(TARGET): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LIBS) -o $(TARGET)

obj/main.o: src/main.cpp
	mkdir -p obj
	$(CC) $(CFLAGS) -c src/main.cpp -o obj/main.o

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
	rm -f $(TARGET).exe

