#ifndef _TIME_MANAGER_H_
#define _TIME_MANAGER_H_

#include "types.h"
#include "Logger.h"
#include <omp.h>
#include <chrono>
#include <unordered_map>
#include <iostream>

/**
 * @brief Clase para controlar el tiempo de la simulacion
 * 
 */
class TimeManager
{
    private:

        Real fps;
        Real lastTimeSave;

        Real time;
        Real startTime;
        Real endTime;
        Real ts;

        std::unordered_map<std::string, Real> intervals;

        Real minTimeStep;
        Real maxTimeStep;

    public:

        TimeManager()
        {
            time = 0.0;
            ts = 0.0001;
        }

        void setTime(Real time) { this -> time = time; }
        void setStartTime(Real startTime) { this -> startTime = startTime; }
        void setEndTime(Real endTime) { this -> endTime = endTime; }
        void setTimeStep(Real ts) { this -> ts = ts; }
        void setFPS(Real fps) { this -> fps = fps; }
        void setMinTimeStep(Real minTs) { minTimeStep = minTs; }
        void setMaxTimeStep(Real maxTs) { maxTimeStep = maxTs; }

        Real getTime() { return time; }
        Real getStartTime() { return startTime; }
        Real getEndTime() { return endTime; }
        Real getTimeStep() { return ts; }
        Real getMinTimeStep() { return minTimeStep; }
        Real getMaxTimeStep() { return maxTimeStep; }

        /*void printOutput()
        {
            struct Interval
            {
                Real start;
                Real stop;
            };

            std::vector<std::string> index;
            unsigned int max;

            std::string label = "Time";
            std::cout << label + std::string(max - label.length(), ' ') << " -> " << time << " s" << std::endl;
            for (const std::string & i: index)
            {
                label = i + std::string(max - i.length(), ' ');
                std::cout << label << " -> " << intervals[i].stop - intervals[i].start << " s" << std::endl;
            }

            index.clear();
            max = 0;
        }*/

        void startCounting(std::string name)
        {
            Real start;

            #ifdef _OPENMP
            start = omp_get_wtime();
            #else
            start = clock() / (Real) CLOCKS_PER_SEC;
            #endif

            intervals[name] = start; 
        }

        void stopCounting(std::string name)
        {
            Real stop;

            #ifdef _OPENMP
            stop = omp_get_wtime();
            #else
            stop = clock() / (Real) CLOCKS_PER_SEC;
            #endif

            LOG(name, " -> ", stop - intervals[name], " s");
        }

        bool hasToSave()
        {
            if (time - lastTimeSave >= 1.0 / fps)
            {
                lastTimeSave = time; // Cambiar esto si luego no coincide el frame final con el que deberia ser respecto al tiempo final
                
                return true;
            }
            else
                return false;
        }

        void clean()
        {
            intervals.clear();
        }
};

#endif
 
