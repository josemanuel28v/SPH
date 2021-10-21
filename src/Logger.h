#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <sstream>
#include <iostream>
//#include <fstream>

class Logger
{
    private:

        // Añadir niveles de log para seleccionar a partir de que nivel se muestran cosas
        // ERROR, WARN, INFO, TIME (importantes info y time para poder ver o no los tiempos)

        static std::ostringstream stream;
        //static std::ofstream stream;

    public:

        template<typename T, typename... Types>
        static void write(T s, Types... args)
        {
            //stream << s;      
            std::cout << s;      

            write(args...);      
        }

        static void write()
        {
            //stream << std::endl;
            std::cout << std::endl;
        }

        static std::string str()
        {
            return stream.str();
        }
};

#define LOG Logger::write

#endif
