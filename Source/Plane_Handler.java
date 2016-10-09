
//----------------------------------Need Working 2D Plane Cases------------------------------

//To do: Plane to Plane

public class Plane_Handler {

			//Calculates Collision point between Sphere and Plane - Moving or Stationary Sphere. Stationary Plane
			
			//Inputs:   Sphere 1 first and final position as array respectively (3 inputs: x,y,z) - P(t)
			//			Start position of Plane centre 3D double (x,y,z)
			//			end position of plane centre 3D double (x,y,z)
			//			Bottom Right Edge of Plane 3D double (x,y,z)
			//			Line 1 from edge reference to start position of edge 3D double (x,y,z)
			//			Line 2 from edge reference to start position of edge 3D double (x,y,z)
			//			Sphere 1 radius
			//
			//Return:	double[3] = (-1,0,0) if no collision
			//
			//Return:	If Collision		-			Point of collision (x,y,z) & time of collision (t) as double[4] respectively
			
			static Collision_Information Linear_Sphere_to_3D_Finite_Plane(double[] Sphere_Initial, double[] Sphere_Velocity, double[][] Plane_Initial, double[] Plane_Velocity, double r)
			{
				
				/*
				 * 		Creates a 2D Box
				 * 		Affine homogeneous transformation Matrix is use to determine required angle of rotation and translate to 
				 * 		get 2D Box parallel to Plane and overlapping Plane.
				 * 		Uses Infinite Plane method to calculate a collision point with infinite plane.
				 * 		If no collision with infinite plane, then collision with finite plane impossible
				 * 		use 2D Box as Plane inputs to generate infinite plane normals.
				 * 		(Plane normal x edges of the Actual plane) = normals perpendicular to plane edges
				 * 		pointing towards the centre of the plane.
				 * 		D = Normal[x*y*z]		P(t) = P(t) - D.N				=>			Reference to origin
				 * 		P(t).N	=	Northo
				 * 		Use Beck's clipping algorithm to determine collision occurs with finite plane.
				 * 		for all orthogonal normals,			Northo.P(t) >= 0		if Northo pointing towards CoM of Plane
				 * 		
				 */
				
				//After testing, cross product of any two lines on a plane has an orthogonal normal
				//to that plane. Thus we can just use the Normal_to_Plane (modded) method to return normal.
				//Will add comment if regression testing proves otherwise.
				//
				//			=>				Above commented and described method not needed to derive normal
				double[] Normal_to_Plane = Maths.Unit_Vector(Maths.Normal_to_Plane(Plane_Initial));
				
				//Derive Plane_Initial -> Plane_Final Translation Matrix
				//Assume Translation Matrix for all points is identical
				double[] Vp = Plane_Velocity;				
				
				/*See comment above. This method is redundant.
				//Use (2D Box normal . actual planes) = 0 & Linear Transformation Matrix Deriving
				//
				//Find normal -| to 2D plane
				//
				// a.b == 0			where a = edge 0 of plane		b = edge n of plane
				//
				// n = a x b
				
				//double[] Normal_to_Plane = Normal_to_Plane(Box_2D);
				*/
				
				//Check Plane Normal -| to all points on actual plane.
				Boolean Normal_to_Plane_Perpendicular_to_Actual_Plane = true;
				
				//Make all plane points with reference to world frame not local frame
				double D = (Normal_to_Plane[0]*Plane_Initial[0][0]) + (Normal_to_Plane[1]*Plane_Initial[0][1]) + (Normal_to_Plane[2]*Plane_Initial[0][2]);
				double[] N_dot_D = Maths.Vector_Multiply(D,Normal_to_Plane);
				
				//Check condition
				for(int i = 0; i < Plane_Initial.length; i++)
				{
					double[] Pact = Maths.Vector_Subtract(Plane_Initial[i],N_dot_D);
					//Print_Vector("Pact",Pact);
					
					if(Maths.Magnitude_3D(Maths.Dot_Product(Normal_to_Plane, Pact)) != 0)
					{
						Normal_to_Plane_Perpendicular_to_Actual_Plane = false;
						System.out.println("Normal not -| to actual plane");
						Maths.Print_Vector("N", Normal_to_Plane);
						return null;
					}
				}
				
				
				//If Code gets Here Means Transformation to infinite plane is correct.
				
				//Check Sphere_to_2D_Infinite_Plane has detected a collision.
				//Use box to get normal -| to infinite plane
				//If returned no collision [-1,0,0] return false (no collision)
				//else continue
				
				Collision_Information C_data = new Collision_Information(1);
				C_data = Linear_Sphere_to_3D_Infinite_Plane(Sphere_Initial, Sphere_Velocity, Plane_Initial, Plane_Velocity, r);
				
				if(C_data == null)
				{
					System.out.println("No Collision Occurred with Infinite Plane");
					return null;
				}
				
				//If got to here, infinite plane collision occurred
				
				
				//At collision time t, calculate position of plane when collision had occurred (for moving plane scenarios)
				//Multiply Plane_Initial by translation matrix of Vp
				//
				// P(t) = P0 + t*Vp
				//
				// Note: For linear motion. Plane normal is constant.
				//
				double[][] Plane_at_Collision = new double[Plane_Initial.length][Plane_Initial[0].length];
				for(int i = 0; i < Plane_Initial.length; i++)
				{
					Plane_at_Collision[i] = Maths.Vector_Add(Plane_Initial[i],Maths.Vector_Multiply(C_data.time,Vp));				//	Collision_Details[4] = t of collision 
				}
				
				//Use normal -| to 2D plane & edges of actual Plane to calculate Northo
				//
				// n x a = Northo		where n = normal -| to 2D Plane.
				//							  a = edge n of actual plane
				//							  Northo is normal -| to actual plane edge pointing towards CoM
				
				double[][] Finite_Plane_Normals = new double[Plane_Initial.length][3];
				Finite_Plane_Normals = Maths.Calculate_Normal_to_Plane_Edges(Plane_at_Collision,Normal_to_Plane);
	
				//Use Beck's Clipping Algorithm with Northo and P(t) collided with infinite plane
				//To decide whether collision has occurred with finite plane.
				//
				// For all Northo ->		Northo.P(t) >= 0
				//
				// Must Satisfy all conditions to be true
				//
				//Reference for each Northo check is taken at world frame. Body frame CoM taken at plane edge
				
				Boolean Collision_with_Finite_Plane = true;
				for(int i = 0; i < Plane_Initial.length; i++)
				{
					double Dact = (Finite_Plane_Normals[i][0]*Plane_at_Collision[i][0]) + (Finite_Plane_Normals[i][1]*Plane_at_Collision[i][1]) + (Finite_Plane_Normals[i][2]*Plane_at_Collision[i][2]);
					double[] Pact = Maths.Vector_Subtract(C_data.Point[0],Maths.Vector_Multiply(Dact,Finite_Plane_Normals[i])); 
					double[] test = Maths.Dot_Product(Finite_Plane_Normals[i],Pact);
					
					if( ((test[0] < 0) || (test[1] < 0)) || (test[2] < 0) )
					{
						System.out.println("Beck's Clipping Algorithm Failed");
						Collision_with_Finite_Plane = false;
						//return null;
					}
				}
				
				//If any of Beck's conditions failed		=>			No collision with finite plane
				//
				if(Collision_with_Finite_Plane == false)
				{
					return null;
				}
				
				//If got here, Collision occurred!!! (Yaaaay)
				
				return C_data;
			}
			
			//Perform 2D surface collision check
			//
			//Returns all edge Normals on same 2D infinite plane pointing towards P0
			//
			//Returns Class containing noramls and points on normal's plane
			//
			//If return null, P0 starts inside finite plane
			
			static Edge_Normal_Information Linear_Sphere_to_2D_Plane_Collision_Normal_Data(double[] Sphere_Initial, double[] Sphere_Velocity, double[][] Plane_Initial, double[] Plane_Velocity, double r)
			{
				//Find 2D Plane line normal is pointing towards P0.
				//
				// Any two lines cross product gives N to plane.
				// To find Nedge, N.Line = Nedge.
				//
				// Once found, perform Sign check to ensure passed through
				//
				//
								
				//Find Plane Normal
				double[] Plane_Normal = Maths.Unit_Vector(Maths.Normal_to_Plane(Plane_Initial));
				
				//Print_Vector("Nplane",Plane_Normal);
				
				//Find normal to plane edges
				//
				// Note function returns normals towards plane CoM
				double[][] Plane_Edge_Normals = Maths.Calculate_Normal_to_Plane_Edges(Plane_Initial,Plane_Normal);
				
				// Find edge's closest to P0 (normals pointing towards P0).
				// If more than 1 edge, below method must be true for at least one of these edges
				//
				// Returns array size = No. Edge Normals with each element having a boolean representing whether the edge normal -> P0
				//
				Boolean[] Edge_Normals_Towards_P0 = Maths.Edge_Normals_Pointing_Towards_P0(Plane_Initial, Plane_Edge_Normals, Sphere_Initial);
				
				//Count number of edge normals -> P0
				int No_Normals_towards_P0 = 0;
				for(int i = 0; i < Edge_Normals_Towards_P0.length; i++)
				{
					if(Edge_Normals_Towards_P0[i] == true)
						No_Normals_towards_P0++;
				}
				
				Edge_Normal_Information Edge_Data = new Edge_Normal_Information(No_Normals_towards_P0);
				
				//If no normals -> P0
				//
				// => P0 Starts collided with finite 2D plane
				if(No_Normals_towards_P0 == 0)
				{
					return null;
				}
				
				// Store Normal Data for N -> P0
				int j = 0;
				for(int i = 0; i < Edge_Normals_Towards_P0.length; i++)
				{
					if(Edge_Normals_Towards_P0[i] == true)
					{
						Edge_Data.Normal_towards_P0[j] = Plane_Edge_Normals[i];
						Edge_Data.Point_on_Normals_Plane[j] = Plane_Initial[i];
						j++;
					}
				}
				
				return Edge_Data;
			}

			//Performs collision check for sphere against infinite 3D Plane using plane radius offset (N.D)
			
			static Collision_Information Linear_Sphere_to_3D_Infinite_Plane(double[] Sphere_Initial, double[] Sphere_Velocity, double[][] Plane_Initial, double[] Plane_Velocity, double r)
			{
				
				//--------------------D tested, Normal direction tested, 
				
				
				// Plane L = {N,D}			=>			N.P + D = 0
				//
				// N = Normal Unit Vector to Plane
				// D = Distance between // plane passing origin and plane L
				// P = Point on Plane
				// r = radius of sphere involved in collision
				//
				// Assuming t= 0 no collision.
				// Sphere collides with plane L with  sphere's surface. Centre of mass of Sphere is aligned on plane L' at collision point.
				//
				// L' = {N, (D-r)}		=>				N.P + D = r
				//
				// plane L' is parallel to plane L at a distance of r away in the direction of N.
				//
				// If P on plane is taken as the point of intercept of sphere and plane L', P(t) = P0 + Vt
				// V = P1 - P0
				//
				// Sub P(t) into L' condition (Stationary Plane)		=>		N.P0 + (Vt).N + D - r = 0
				// Rearrange to solve for t:									t = | (r - D - N.P0) / V.N |
				//
				// For moving planes									=>		D(t) = D0 + Vt
				// Here:														Vd = D1 - D0
				//																N.P(t) + D(t) - r = 0
				//																N.P0 + (Vt).N + D0 + Vd.t - r = 0
				//																t(V.N + Vd) = r - D0 - N.P0
				// Rearrange to solve for t:									t = | (r - |D0 - N.P0| ) / ( |Vp.N + Vd| ) |
				//
				// Note D0 taken with reference to plane
				//
				// As P(t_collision) is collision point between Sphere CoM and L'. Collision point of sphere surface and L 
				// is +N.r away from P(t_collision) as
				//
				// L' = {N, D - r}			&&			L = {N, D}				=> -r difference thus +r needed	
				//
				// => C = P(t) + N.r
				//
				// Where C = first point of collision with plane L and Sphere.
				//
				// Non-Collision Conditions:			(N.V = 0) &&  (N.P0 + D = 0)		- V parallel to Normal and is not penetrating plane
				//										! 0 <= t <= 1
				
				// Find normal N from plane lines
				// calculate:		N		L1		L2		D		N.P0		N.V			N.r			N.V
				double[] edge_start = Plane_Initial[0];
				double[] edge_end = Maths.Vector_Add(Plane_Initial[0],Plane_Velocity);
				
				double[] Vp = Sphere_Velocity;
				
				System.out.println("----");
				double[] N = Maths.Unit_Vector(Maths.Normal_to_Plane(Plane_Initial));


				// Check normal is pointing in direction of sphere's starting point (assuming ds -> 0 for L)
				// Normalised normal is taken as vector with origin as reference point (world translated by -Q)
				// To normalise point with respect to plane: Simply minus plane reference pt (edge) from point.
				//
				// D = Normal*[x,y,z]
				// P_act = P(t) - N.D;
				//
				// P_act = point position w. ref. to plane
				// D = distance from plane to // plane passing through origin
				//
				// condition:		N_dot_Pact >= 0				-			plane (reference) to point vector is in same direction to normal N
				//												else		plane (reference) to point vector is in opposite direction to normal N
				//
				//					Pact = point P with reference to plane edge where plane edge is world origin
				//					N = normalised normal of plane.
				//					Q = plane reference point (edge in calculation)
				//
				// N Must be in direction of point (reference) to plane (N pointing away from point) to calculate correct t - same direction as Vd
				// 
				// Proof:				C = P(t) - r.N				-				
				//
				// P(0) is on same side of plane as P(1) thus P(t). If P(t) is on left side of plane L, theoretical plane, L' is on left side of plane L.
				// delta between L' and L is r.N where N is in the direction of P(t) pointing towards plane L [normal pointing away from P(0)].
				// In this scenario, N is pointing to the right hand side of the plane.
				//
				// Condition must stand for all components of N_dot_Pact
				//
				
				// D = Normal*[x,y,z]
				// P_act = P(t) - N.D;
				//
				// P_act = point position w. ref. to plane
				// D = distance from plane to // plane passing through origin
				//
				double D = (N[0]*edge_start[0]) + (N[1]*edge_start[1]) + (N[2]*edge_start[2]);
				double[] N_dot_D = Maths.Vector_Multiply(D,N);
				
				double[] Pact = Maths.Vector_Subtract(Sphere_Initial, N_dot_D); 
				double[] N_dot_Pact = Maths.Dot_Product(N,Pact);
				
				if( ( (N_dot_Pact[0] <= 0) && (N_dot_Pact[1] <= 0) ) && (N_dot_Pact[2] <= 0) )				// If N is pointing away from P, reverse N's direction
				{
					//System.out.println("Invert N");
					N = Maths.Vector_Multiply(-1,N);
				}

				//Preliminary Step:				N.(P(0)/|P(0)|) == N.(P(1) / |P(1)|) for each plane component		=>			No Collision
				// Couldn't Get this code to work.
				// Tbh I think the actual computation (taking correct reference points, doing basic algebra) is not faster than N.V and 0<t<1 methods presented
				
				
				double[] N_dot_P0 = Maths.Dot_Product(N, Sphere_Initial);
				double[] N_dot_Vp = Maths.Dot_Product(N, Vp);
				double[] N_dot_r = Maths.Vector_Multiply(r,N);
				
				// Calculate D
				// D = Ax + By + Cz				-			Equation of a plane where N = [A, B, C] of a normalized vector
				//
				// As parallel plane to L passes through origin, plane O;		D = Nx(0) + Ny(0) + Nz(0)
				//																D = 0
				// Distance from plane passing through origin to origin is 0.
				// Both planes have same right hand side of the equation
				// As plane L is // to plane 0. Has same right normal thus same right hand side of the equation.
				//
				// for L:					D = Nx(X) + Ny(Y) + Nz(Z)		where		X,Y,Z are co-ordinates of any point on plane L eg. edge
				//
				// Using edge as X,Y,Z.	Normal for moving plane is constant. 
				// D0 is calculated using Edge at t = 0 (edge_start)
				// D1 is calculated using Edge at t = 1 (edge_end)
				//
				// Vector Vd describes plane's linear motion.
				// Vd = P(1) - P(0)							-				P is taking at edge's point in space.
				//
				// Decide D0 sign depending on origin's position with reference to plane
				//
				// Vd is with reference to the end point
				//
				
				
				double D0 = N[0]*edge_start[0] + N[1]*edge_start[1] + N[2]*edge_start[2]; 
				double D1 = N[0]*edge_end[0] + N[1]*edge_end[1] + N[2]*edge_end[2]; 
				
				//As D is calculated initially as magnitude find correct sign.
				//I don't know why this condition stands, but the proof work for a and 1 for t=0.5 collision at back of journal proves.

			//	Vd = Vector_Subtract(edge_end,edge_start);
				double Vd = D1 - D0;
				
				//First Condition Check:		N.V = 0		&&		N.P0 + D = 0		->		Parallel V and no intersect with plane L
				//
				// Point Vector & plane may be in the same infinite plane => Either 2D occurrance
				//
				// If not in the same plane, then implies no collision
				double[][] Two_D_Scenario_Normals = null;
				Edge_Normal_Information M = null;
				
				Maths.Print_Vector("Vp",Vp);
				Maths.Print_Vector("N",N);
				
				double[] Sphere_Final = Maths.Vector_Add(Sphere_Initial, Sphere_Velocity);
				double[][] Plane_Final = new double[Plane_Initial.length][3];
				
				int i;
				for(i = 0; i< Plane_Initial.length; i++)
				{
					Plane_Final[i] = Maths.Vector_Add(Plane_Initial[i],Plane_Velocity);
				}
				
				if(	(Maths.Magnitude_3D(N_dot_Vp) == 0) &&  (  (Maths.Magnitude_3D(N_dot_P0) + D0) == 0 ) )
				{
					//Perform 2D Check
					//
					// IF plane is 2D, then 1 dimension must equal to 0 for Plane & Vp. 
					// If this is the case call 2D plane collision method
					//
					// if plane is not 2D then return null.
					
					//Check 2D plane
					//
					// If 2D: Normal to plane, N, is also the normal to Vp.
					//
					// Dp0 = | N*[P0] |
					// Dq0 = | N*[Plane0] |
					// Dp1 = | N*[P1] |
					// Dq1 = | N*[Plane1] |
					// Dqn == Dpn for all n			=>			2D Collision Scenario
					
					double Dp0 = Maths.Magnitude_3D(Maths.D_Plane_to_Parallel_Plane_Through_Origin(N,Sphere_Initial));
					double Dq0 = Maths.Magnitude_3D(Maths.D_Plane_to_Parallel_Plane_Through_Origin(N,Plane_Initial[0]));
					double Dp1 = Maths.Magnitude_3D(Maths.D_Plane_to_Parallel_Plane_Through_Origin(N,Sphere_Final));
					double Dq1 = Maths.Magnitude_3D(Maths.D_Plane_to_Parallel_Plane_Through_Origin(N,Plane_Final[0]));

					System.out.println("Dp0: " + Dp0);
					System.out.println("Dq1: " + Dp1);
					System.out.println("Dq0: " + Dq0);
					System.out.println("Dq1: " + Dq1);
					
					// If distance between both plane positions to origin and both point positions to origin is same
					// => 2D Scenario (Point, Vp and Plane, are on same infinite plane!!!!!!! i.e. [x,y] share same z)
					if( (Dp0 == Dp1)  && (Dp0 == Dq0) && (Dp0 == Dq1))
					{
						System.out.println("2D Collision Scenario");
						
						//M holds information about the 2D Normals Pointing towards P0
						M = Linear_Sphere_to_2D_Plane_Collision_Normal_Data(Sphere_Initial,Sphere_Velocity,Plane_Initial,Plane_Velocity, r);
						
						if(M.Normal_towards_P0.length != 0)
							Two_D_Scenario_Normals = new double[M.Normal_towards_P0.length][3];
					}
					else
					{
						//Point Starts inside finite plane
						//If No 2D Collision then no collision occurs
						System.out.println("No Collision (N.Vp = 0): Plane and Vector are parallel");
						return null;
					}
				}
				
				/*
				 * If Not // to plane must pass through plane at some point. 
				 * 
				 * It is actually quicker to calculate t @ collision and comapre against 0 <= t <= 1 than it is
				 * to calculate whether sign changes (point passes through plane) then proceed with t @ collision
				 * calculation anyway. This is because t @ collision only needs 1 if statement (check) to see
				 * whether it has passed through within that collision frame. If 0 <= t <= 1, then collision is 
				 * not within animation frame!!!
				 */
				
				//Calculate t of collision point with @L'

				
				
				//double t = Magnitude_Scalar(((r -  D0) - Magnitude_3D(N_dot_P0) )/ (Magnitude_3D(N_dot_V) + V));
				double t_collision = 1;
				
				if(Two_D_Scenario_Normals == null)
					t_collision = Maths.Magnitude_Scalar(((r +  D0) - Maths.Magnitude_3D(N_dot_P0) )/ (Maths.Magnitude_3D(N_dot_Vp) + Vd));
				else
				{
					//As N for 3D is taken in direction of point, N 2D is calculated towards CoM. Invert direction as know -1.N will point towards P0
					
					//Will check time of collision against all normals pointing towards P0
					//Lowest time in collision will be collision time
					for(i = 0; i < M.Normal_towards_P0.length; i++)
					{
						N = M.Normal_towards_P0[i];
						N = Maths.Vector_Multiply(-1,N);
						D0 = Maths.Magnitude_3D(Maths.D_Plane_to_Parallel_Plane_Through_Origin(N,M.Point_on_Normals_Plane[i]));
					
						double t = Maths.Magnitude_Scalar(((r +  D0) - Maths.Magnitude_3D(N_dot_P0) )/ (Maths.Magnitude_3D(N_dot_Vp) + Vd));
						
						if(t < t_collision && t >= 0)
							t_collision = t;
					}
				}
				
				System.out.println("t: " + t_collision);
				
				//Second Collision Check:		Collision occurs in animation interval
				if( (t_collision < 0) || (t_collision > 1) )
				{
					System.out.println("No Collision !(0<t<1)");
					return null;
				}
					
				//---------------------If its made it this far, collision in animation frame does occurs-----------------------------------
				
				// Calculate point of collision, P(t) with L'
				// P(t_collision) = t.V + P0			-		V = P1 - P0
				double[] P = Maths.Vector_Add(Maths.Vector_Multiply(t_collision,Vp), Sphere_Initial);
				
				// Calculate point of collision, C(t) with L
				// C(t_collision) = P(t_collision) - N.r
				// 
				double[] C = Maths.Vector_Subtract(P, N_dot_r);
				
				Maths.Print_Vector("P along L': ", P);
				Maths.Print_Vector("C along L: ", C);	
								
				Collision_Information C_data = new Collision_Information(1);
				C_data.Point[0] = C; 
				C_data.time = t_collision;
				
				return C_data;
			}
			
			static Collision_Information Linear_Point_to_3D_Plane(double[] Sphere_Initial, double[] Sphere_Velocity, double[][] Plane_Initial, double[] Plane_Velocity)
			{
				return Linear_Sphere_to_3D_Finite_Plane(Sphere_Initial, Sphere_Velocity, Plane_Initial, Plane_Velocity, 0);
			}
			
			static Collision_Information Linear_1D_Vector_to_3D_Plane(double[] Sphere_Initial, double[] Sphere_Velocity, double[][] Plane_Initial, double[] Plane_Velocity)
			{
				return Linear_Sphere_to_3D_Finite_Plane(Sphere_Initial, Sphere_Velocity, Plane_Initial, Plane_Velocity, 0);
			}
}
