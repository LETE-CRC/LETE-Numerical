# MovMesh
This is a simple Python routine to perform the 3D deformation of a file generated using OpenFOAM.

To run it, it is necessary either have an OpenFOAM case prepared or download the full folder from Github and, then, simples type in a terminal console "python basicDiffusion.py".

Results will be written at folder "1" and can be seen on ParaView by typing in a terminal console "paraFoam".

The default mesh is a cube with 0.1m of size and with 125 volumes. The boundaries are always positionated in such way that one of the axis x, y and z is normal to it. TOP WALL and BOTTOM WALL have the y axis as their normal, but the top wall is at y+ and the bottom wall is at y-. Following this logic, EAST WALL is on x+ while WEST WALL is on x- and NORTH WALL and SOUTH WALL are, respectivily, on z+ and z-. Some of those boundaries are show on the skecth bellow.

					       * * * * * * * * * * * *
					      *                     **
					     *                     * *
					    *        TOP          *  *
					   *         WALL        *   *
					  *                     *    *
					 *                     *     *
					* * * * * * * * * * * * EAST *
					*                     * WALL *
					*                     *      *
					*                     *      *
					*                     *     *
					*     NORTH WALL      *    *
					*                     *   *
					*                     *  *
					*                     * *
					*                     **
					* * * * * * * * * * * *
