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

        struct Interval
        {
            Real start;
            Real stop;
        };

        Real fps;
        Real lastTimeSave;

        unsigned int frame;
        Real time;
        Real startTime;
        Real endTime;
        Real ts;

        std::unordered_map<std::string, Interval> intervals;

        Real minTimeStep;
        Real maxTimeStep;

    public:

        TimeManager()
        {
            time = 0.0;
            ts = 0.0001;
            frame = 1;
        }

        void setFrame(unsigned int frame) { this -> frame = frame; }
        void setTime(Real time) { this -> time = time; }
        void setStartTime(Real startTime) { this -> startTime = startTime; }
        void setEndTime(Real endTime) { this -> endTime = endTime; }
        void setTimeStep(Real ts) { this -> ts = ts; }
        void setFPS(Real fps) { this -> fps = fps; } // Se debe llamar antes que setMinTimeStep()
        void setMinTimeStep(Real minTs) { Real minFPSTimeStep = 1.0 / fps; minFPSTimeStep < minTs? minTimeStep = minFPSTimeStep : minTimeStep = minTs; }
        void setMaxTimeStep(Real maxTs) { maxTimeStep = maxTs; }

        Real getFrame() { return frame; }
        Real getTime() { return time; }
        Real getStartTime() { return startTime; }
        Real getEndTime() { return endTime; }
        Real getTimeStep() { return ts; }
        Real getFPS() { return fps; }
        Real getMinTimeStep() { return minTimeStep; }
        Real getMaxTimeStep() { return maxTimeStep; }


        void startCounting(std::string name)
        {
            Real start;

            #ifdef _OPENMP
            start = omp_get_wtime();
            #else
            start = clock() / (Real) CLOCKS_PER_SEC;
            #endif

            intervals[name].start = start; 
        }

        void stopCounting(std::string name)
        {
            Real stop;

            #ifdef _OPENMP
            stop = omp_get_wtime();
            #else
            stop = clock() / (Real) CLOCKS_PER_SEC;
            #endif

            intervals[name].stop = stop; 

            LOG(name, " -> ", stop - intervals[name].start, " s");
        }

        Real getInterval(std::string name)
        {
            return intervals[name].stop - intervals[name].start;
        }

        bool hasToSave()
        {
            if (time - lastTimeSave >= 1.0 / fps)
            {
                lastTimeSave = (Real) frame / fps;
                
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
 
