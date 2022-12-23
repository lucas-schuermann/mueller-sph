#######################################################################################################

# Mac OS X
SPH_INCLUDE_PATH      = -I/usr/local/include/
SPH_LIBRARY_PATH      = -L/usr/local/lib/
SPH_OPENGL_LIBS       = -framework OpenGL -framework GLUT

# # Linux
# SPH_INCLUDE_PATH      =
# SPH_LIBRARY_PATH      =
# SPH_OPENGL_LIBS       = -lglut -lGL -lX11

# # Windows / Cygwin
# SPH_INCLUDE_PATH      = -I/usr/include/opengl
# SPH_LIBRARY_PATH      = -L/usr/lib/w32api
# SPH_OPENGL_LIBS       = -lglut32 -lopengl32

#######################################################################################################

TARGET = sph
CC = g++
LD = g++
CFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated -pedantic -Wno-vla-extension $(SPH_INCLUDE_PATH) -I./include -I./src -DNDEBUG
LFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated -Werror -pedantic $(SPH_LIBRARY_PATH) -DNDEBUG
LIBS = $(SPH_OPENGL_LIBS)

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
