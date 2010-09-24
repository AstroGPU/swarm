/*************************************************************************
 * Copyright (C) 2010 by Mario Juric  and the Swarm-NG Development Team  *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 3 of the License.        *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the                         *
 * Free Software Foundation, Inc.,                                       *
 * 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ************************************************************************/

/*! \file cpu_ensemble.h
 *   \brief Declares cpu_ensemble class 
 *
*/

/// The main namespace for the Swarm-NG library
namespace swarm {
/**
        \brief Specialization of ensemble, with data stored on the host (CPU)

        Use this class to create/read/write ensembles in host memory.
*/
class cpu_ensemble : public ensemble
{
        public:
		//constructors
                cpu_ensemble();                                 //!< instantiate an empty ensemble
                cpu_ensemble(int nsys, int nbod);               //!< instantiate an ensemble with room for nsys systems with nbod each
                cpu_ensemble(const gpu_ensemble &source);       //!< instantiate a copy of source ensemble
                cpu_ensemble(const cpu_ensemble &source);       //!< instantiate a copy of source ensemble

                ~cpu_ensemble() { free(); }

                //Methods defined in cpu_ensemble.cpp
                void copy_from(const cpu_ensemble &source);     //!< Copy the data from the CPU
                void copy_from(const gpu_ensemble &source);     //!< Copy the data from the GPU

                void free();                                            //!< Deallocate CPU memory
                int pack();
                void replace_inactive_from(cpu_ensemble &src, const int offset);
                void reset(int nsys, int nbod, bool reinitIndices = true);      //!< Allocate CPU memory for nsys systems of nbod planets each


		//Inline methods
                void set_active(const std::vector<int>& keep_flag)
                {
                  for(int sysid=0;sysid<nsys();++sysid)
                    if(keep_flag[sysid]) ensemble::set_active(sysid);
                };

                void set_inactive(const std::vector<int>& halt_flag)
                {
                  for(int sysid=0;sysid<nsys();++sysid)
                    if(halt_flag[sysid]) ensemble::set_inactive(sysid);
                };

                void set_active(const int sys) { ensemble::set_active(sys); }
                void set_inactive(const int sys) { ensemble::set_inactive(sys); }


        private:
                cpu_ensemble &operator=(const cpu_ensemble&);   //!< disallow copying
};
//end class cpu_ensemble
}
