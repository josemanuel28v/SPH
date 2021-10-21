#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/**
 * Class InputParser 
 */
class InputParser
{
    public:

        /**
         * @brief Inicializa un objeto de la clase a partir del array de argumentos y el numero de argumentos
         * 
         * @param argc Numero de argumentos
         * @param argv Array de argumentos
         * 
         * Almacena los argumentos del array argv[] en un vector de strings
         */
        InputParser (int argc, char *argv[])
        {
            for (int i = 1; i < argc; ++i)
                tokens.push_back(string(argv[i]));
        }

        /**
         * @brief Devuelve el argumento que acompaña a la opcion recibida como parametro
         * 
         * @param option Opcion para la cual se busca el argumento
         * 
         * @return El argumento buscado o una string vacia si option no se encuentra en el vector de tokens
         * 
         * Busca la string option en el vector de tokens y si lo encuentra devuelve el token siguiente
         */ 
        static string getCmdOption(string option) 
        {
            vector<string>::const_iterator itr;
            itr =  find(tokens.begin(), tokens.end(), option);

            if (itr != tokens.end() && ++itr != tokens.end())
                return *itr;
            else
                return string();
        }

        /**
         * @brief Devuelve true si se encuentra la string option en el vector de tokens
         * 
         * @retval True si se encuentra el token, false si no se encuentra
         */
        static bool cmdOptionExists(string option) 
        {
            return find(tokens.begin(), tokens.end(), option) != tokens.end();
        }

        /**
         * @brief Devuelve la cantidad de tokens
         * 
         * @return La cantidad de tokens
         */
        static unsigned size()
        {
            return tokens.size();
        }

    private:

        static vector<string> tokens;      ///< Vector con las cadenas del array de argumentos
};