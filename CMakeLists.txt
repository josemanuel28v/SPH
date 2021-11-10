PROJECT(SPH)
CMAKE_MINIMUM_REQUIRED(VERSION 3.11) # Indicamos la versión mínima que necesitamos de CMake

# Incluir directorios con codigo a compilar
include_directories(${SPH_SOURCE_DIR}/extern ${SPH_SOURCE_DIR}/src)

# Establecer los flags para la compilacion
SET(CMAKE_CXX_FLAGS "-O3 -mavx -ffast-math -std=c++17 -march=native -mtune=native -Wall -Wpedantic -fopenmp -lglut -lGL -lGLEW")

# Todos los cpp para los que se generara un .o
FILE(GLOB SPH_SRC "${SPH_SOURCE_DIR}/Extern/*.cpp"
                  "${SPH_SOURCE_DIR}/src/*.cpp")

ADD_EXECUTABLE(SPH ${SPH_SRC}) # Crea un ejecutable llamado SPH en Linux o SPH.exe en Windows.

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)

include_directories(GL glut GLEW)

target_link_libraries(SPH GL glut GLEW)
