#ifndef _TIME_MANAGER_H_
#define _TIME_MANAGER_H_

#include "types.h"
#include <omp.h>
#include <chrono>

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
        Real ts;

        Real start;
        Real stop;

        Real minTimeStep;
        Real maxTimeStep;

    public:

        TimeManager()
        {
            time = 0.0;
            ts = 0.0001;

            start = 0.0;
            stop = 0.0;
        }

        void setTime(Real time) { this -> time = time; }
        void setTimeStep(Real ts) { this -> ts = ts; }
        void setFPS(Real fps) { this -> fps = fps; }
        void setMinTimeStep(Real minTs) { minTimeStep = minTs; }
        void setMaxTimeStep(Real maxTs) { maxTimeStep = maxTs; }

        Real getTime() { return time; }
        Real getTimeStep() { return ts; }
        Real getInterval() { return stop - start; }
        Real getMinTimeStep() { return minTimeStep; }
        Real getMaxTimeStep() { return maxTimeStep; }

        void startCounting()
        {
            #ifdef _OPENMP
            start = omp_get_wtime();
            #else
            start = clock() / CLOCKS_PER_SEC;
            #endif
        }

        void stopCounting()
        {
            #ifdef _OPENMP
            stop = omp_get_wtime();
            #else
            stop = clock() / CLOCKS_PER_SEC;
            #endif
        }

        bool hasToSave()
        {
            if (time - lastTimeSave >= 1.0 / fps)
            {
                lastTimeSave = time;
                
                return true;
            }
            else
                return false;
        }
};

#endif
 
