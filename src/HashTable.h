#ifndef _HASHTABLE_H_
#define _HASHTABLE_H_

#include <vector>
#include <unordered_map>
#include <glm/glm.hpp>
#include <algorithm>

#include "types.h"

using std::vector;
using std::unordered_map;

/**
 * Estructura que contiene la funcion de hash utilizada
 * en la estructura unordered_map 
 */
struct Hash
{
    /**
     * @brief Calcula el hash del parametro recibido
     * 
     * @param v Coordenadas discretizadas que representan la celda de la particula
     * @return El hash correspondiente al parametro v recibido
     */
    size_t operator()(Vector3i v) const
    {
        return v.x * 73856093 + v.y * 19349663 + v.z * 83492791;
    }
};

/*
 *  Class Hash Table
 */
class HashTable
{
    public:

        struct pointInfo
        {
            unsigned int id;
            unsigned int pointSetId;
        };

        struct Bucket                                   ///< Estructura que contiene las particulas y los vecinos de cada celda. Es el value de la pareja {key, value}
        {                                               
            vector<pointInfo> points;                     
            vector<pointInfo> neighbors;                    
            vector<pointInfo> boundary_points;            
            vector<pointInfo> boundary_neighbors;           
        };   
        /*
        En lugar de que el struct Bucket sea el tipo de dato del valor, que sea un unsigned que indique el indice de celda??  y tener las celdas en un vector aparte
        */                                           
                                                        
        unordered_map<Vector3i, Bucket, Hash> table;       ///< Tabla hash que utiliza como clave una celda y como valor el bucket (que contiene las particulas y los vecinos de esa celda) 
        vector<Vector3i> keys;

        Real h;                                       ///< Distancia de suavizado (Smoothing length)

        HashTable() {}

        /**
         * @brief Inicializa la tabla hash reservando el espacio correspondiente
         * 
         * @param size El espacio que se va a reservar
         */
        HashTable(unsigned size)
        {
            reserve(size);
        }

        /**
         * @brief Inserta el id de una particula en la tabla hash
         * 
         * @param pos Posicion de la particula
         * @param id Id de la particula
         * 
         * El id de la particula se inserta en el vector points del bucket correspondiente a la celda obtenida de la discretizacion de su posicion
         * La posicion se discretiza dividiendo cada componente de la posicion entre la distancia de suavizado
         */
        void insertFluidParticle(Vector3r pos, unsigned int id, unsigned int pointSetId)
        {
            Vector3i key = floor(pos / h);

            pointInfo value = {id, pointSetId};
            //unsigned value = id;

            table[key].points.push_back(value);
        }

        /**
         * @brief Inserta cada boundary particle recibida en el vector boundary_points de su bucket correspondiente
         *  
         * La insercion de las boundary particles solo se realiza una vez al inicio de la simulacion
         */
        void insertBoundaryParticle(Vector3r pos, unsigned id, unsigned int pointSetId)
        {
            Vector3i key = floor(pos / h);

            pointInfo value = {id, pointSetId};
            //unsigned value = id;

            table[key].boundary_points.push_back(value);
        }

        /**
         * @brief Encuentra las particulas vecinas correspondientes a cada celda
         * 
         * Para cada celda (utilizando la clave de una particula vecina prototipo (neight_key) se exploran las celdas de su alrededor (incluida la propia)
         * Los identificadores de las particulas que se encuentran en estas celdas vecinas se almacenan en el vector neights del bucket correspondiente a la
         * celda para la cual explorabamos sus celdas vecinas
         */
        #ifndef _OPENMP
        void neighborhoodSearch()
        {
            for (auto& cell: table)
                for (int x = - 1; x <= 1; x++)
                    for (int y = - 1; y <= 1; y++)
                        for (int z = - 1; z <= 1; z++)
                        {
                            Vector3i dir(x, y, z);
                            Vector3i neigh_key = cell.first + dir;
                            auto neigh_cell = table.find(neigh_key);
                            
                            if (neigh_cell != table.end())
                            {
                                cell.second.neighbors.insert(cell.second.neighbors.end(), 
                                                          neigh_cell -> second.points.begin(), 
                                                          neigh_cell -> second.points.end()); 

                                cell.second.boundary_neighbors.insert(cell.second.boundary_neighbors.end(),       // De esta manera se consigue tener las boundary particles que son vecinas de una celda en un vector distinto de neighbors
                                                          neigh_cell -> second.boundary_points.begin(), 
                                                          neigh_cell -> second.boundary_points.end()); 
                            }
                        }
        }
        #else
        void neighborhoodSearch()
        {
            // Asi mejora respecto a arriba pero oes necesario utilizar como tabla hash un array donde cada elemento sea la celda en la que esta cada particula
            vector<Vector3i> keys;
            for (auto& cell: table)
                keys.push_back(cell.first);

            #pragma omp parallel for
            for (auto key: keys)
            {
                auto & cell = table[key];
                for (int x = - 1; x <= 1; x++)
                    for (int y = - 1; y <= 1; y++)
                        for (int z = - 1; z <= 1; z++)
                        {
                            Vector3i dir(x, y, z);
                            Vector3i neigh_key = key + dir;
                            auto neigh_cell = table.find(neigh_key);
                            
                            if (neigh_cell != table.end())
                            {
                                cell.neighbors.insert(cell.neighbors.end(), 
                                                          neigh_cell -> second.points.begin(), 
                                                          neigh_cell -> second.points.end()); 

                                // Si las partículas son estáticas siempre van a ser vecinas de los mismos vóxeles por lo que esto se podría hacer solo una vez
                                cell.boundary_neighbors.insert(cell.boundary_neighbors.end(),       // De esta manera se consigue tener las boundary particles que son vecinas de una celda en un vector distinto de neighbors
                                                          neigh_cell -> second.boundary_points.begin(), 
                                                          neigh_cell -> second.boundary_points.end()); 
                            }
                        }  
            }   
        }
        #endif

        /**
         * @brief Limpia los vectores points y neighbors de cada bucket en la tabla hash
         * 
         */
        void clear()
        {
            for (auto& cell: table)
            {
                cell.second.points.clear();
                cell.second.neighbors.clear();
                cell.second.boundary_neighbors.clear();
            }
        }

        /**
         * @brief Limpia el vector boundary points de cada bucket en la tabla hash
         * 
         */
        void clearB()
        {
            for (auto& cell: table)
                cell.second.boundary_points.clear();
        }

        /**
         * @brief Reserva espacio en memoria para la tabla hash
         * 
         * @param size Espacio en memoria que se va a reservar
         */
        void reserve(size_t size)
        {
            table.reserve(size);
        }

        /**
         * @brief Devuelve el tamaño de la tabla
         * 
         * @return El tamaño de la tabla
         */
        int size()
        {
            return table.size();
        }

        /**
         * @brief Establece la longitud de suavizado 
         * 
         * @param h Longitud de suavizado
         * 
         * El parametro h sera utilizado para discretizar las posiciones de las particulas que se van a insertar
         */
        void setSmoothingLength(Real h)
        {
            this -> h = h;
        }
};

#endif