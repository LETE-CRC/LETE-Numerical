# MovMesh
This is a simple Python routine to perform the 3D deformation of a file generated using OpenFOAM.

To run it, it is necessary either have an OpenFOAM case prepared or download the full folder from Github and, then, simples type in a terminal console "python MovMesh.py".

Results will be written at folder "1" and can be seen on ParaView by typing in a terminal console "paraFoam".

The default mesh is a cube with 0.1m of size and with 125 volumes. The boundaries are always positionated in such way that one of the axis x, y and z is normal to it. TOP WALL and BOTTOM WALL have the y axis as their normal, but the top wall is at y+ and the bottom wall is at y-. Following this logic, EAST WALL is on x+ while WEST WALL is on x- and NORTH WALL and SOUTH WALL are, respectivily, on z+ and z-. The nodes used to define each boundaries are shown on the skecth bellow and the table of nodes/boundaries correlation is written bellow the figure.

|-|-|-|-|-|-|-|-|-| FIGURE 1 - GEOMETRY SKETCH  |-|-|-|-|-|-|-|-|
|		       * * * * * * * * * * * *N6(0.1, 0.1, 0.1)	|
|		      *.N7(0, 0.1, 0.1)     **			|
|		     * .                   * *			|
|		    *  .                  *  *			|
|		   *   .                 *   *			|
|		  *    .                *    *			|
|		 *     .               *     *			|
|	        * * * *.* * * * * * * *      *			|
|N3(0, 0.1, 0)	*      .            N2*      *			|
|		*      . (0.1, 0.1, 0)*      *			|
|		*      . . . . . . . .*. . . *N5(0.1, 0, 0.1) 	|
|		*     .N4(0, 0, 0.1)  *     *			|
|		*    .                *    *			|
|		*   .                 *   *			|
|		*  .                  *  *			|
|		* .                   * *			|
|		*.                    **			|
|		* * * * * * * * * * * *				|
|           N0(0,0,0)             N1(0.1,0,0) 			|
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

|-|-|TABLE 1 - BOUNDARIES DEFINITION|-|-|
|	BOUNDARY|NODES			|
|	 top	|2 6 7 3		|
|	 bottom	|0 4 5 1		|
|	 north	|7 6 5 4		|
|	 south	|0 1 2 3		|
|	 east	|0 3 7 4		|
|	 west	|1 5 6 2		|
|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
