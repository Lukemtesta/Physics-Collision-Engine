

//--------------------------------------------------Updated Box Methods 24.2.13------------------------------------------


import java.math.BigDecimal;
import java.math.MathContext;

/* Bounding Box (Non-Angular):
	 * 
	 * Box to Box (Effective Radius) (Main)
	 * Box to Sphere Y (Main)
	 * Circumscribed Box to Plane (main)
	 * Circum Box to Plane (main)
	 * Accurate Box to Plane (main)
	 * Box to Point
	 * Box to 1D Plane 
	 * Box to 2D Plane (Scale Rotate/Scale Box into plane. Box to Box)  (Main)
	 */ 

public class Box_Handler {

	
	//Method to calculate the collision between a box and a plane - CIRCUMSCRIBED SPHERE = RADIUS
	//
	//Inputs: Box Start, Box End, Plane Start, Plane End
	//
	//Return:		null no collision
	//				Collision Information Class if collided containing details of collision
	//
	//
	// Note: Accurate for Single Vertex and Edge Collision Only!
	//		 Why? Because face collision is RST not cricumscribed radius apart. At this point offsets plane
	//		by difference which may not be in plane's original vector path
	
	static Collision_Information Linear_Circumscribed_Box_to_Plane(double[][][] Box0, double[] Box_Velocity, double[][] Plane0, double[] Plane_Velocity)
	{
		Collision_Information C = new Collision_Information(1);
		
		/*This method approximates a cube as 8 vertices. 
		 * 
		 * Cube can only collide 3 ways: 		Single Corner
		 * 										Dual corners (edge collision)
		 * 										4 corners (face collision)
		 * 
		 * Effective radius is the distance between the CoM and a point on the edge of a cube
		 * 
		 * Reff is either inscribed sphere (internal sphere with radius Reff connected to mid-point of cube faces)
		 * 		or Circumscribed Sphere (external Sphere intersects with corners of cube)
		 * 
		 * In this case circumscribed sphere radius, Reff is used.
		 * 
		 * Reff = (sqrt(3)/2) * (|R.N| + |S.N| + |T.N|)
		 * 
		 * where [R,S,T] is length of edge where an endpoint is referenced at world frame origin
		 * thus:				[R, S, T] = [(sqrt(3)/2) * (-1) ^ n & 1, (sqrt(3)/2) * (-1) ^ n & 2, (sqrt(3)/2) * (-1) ^ n & 4]
		 * 
		 * and N = normal to edge pointing away from cube
		 * 
		 * or [R,S,T] = CoM - Box[0][0];
		 * 
		 * Infinite Plane L = {N,D} where N.P + D = 0
		 * P = point of object touching plane
		 * r = radius of object touching plane
		 * 
		 * L' is an infinite plane parallel to plane L offset by r. 
		 * L' = {N, D - r}
		 * 
		 * Check CoM passes through L':		NL'.Q(0)/|NL'.Q(0)| != NL'.Q(1)/|NL'.Q(1)|
		 * 
		 * Where NL' = Normal of infinite plane L'
		 * 
		 * 
		 * Find time L' intersects with CoM
		 * 
		 * NL'.Q(t) + D = r
		 * NL'.Q(0) + NL'.Vq.t + D = r
		 * t = [r - NL'.Q(0)]/Vq.NL'
		 * 
		 * Where Q(t) = Q(0) + Vq.t
		 * 		 Vq = Q(1) - Q(0)
		 * 
		 * If 0 <= t <= 1 then collision occurs. Else return null.
		 * 
		 * Zn = CoM(t) + [ +/- R +/- S +/- T]
		 * Where Zn - Vertex n position on collision
		 * 
		 * Check orientation of box on collision. Find signs of [R,S,T] such that NL'.Zn -> -oo
		 * The result of the dot product must be negative. Thus is positive result, invert sign.
		 * Special case is dot product is = to 0. 
		 * 
		 * Once [R,S,T] signs found, do dot product with normal to cube edge to find cube's orientation.
		 * 
		 * No dot product = 0 			=>			Single Vertex Collision
		 * 1 dot product = 0			=>			Edge Collision (Two Vertex)
		 * 2 Dot product = 0			=>			Face Collision ( 4 Vertex)
		 * 3 Dot Product = 0			=>			Plane is cube and aligned with Cube
		 * 
		 * Dot Products are |N.R|, |N.S| & |N.T|.
		 * Where N is normal to cube's edge
		 * 
		 * C = Q(t) - 0.5*[sgnR(R.N) + sgnS(S.N) + sgnT(T.N)]
		 * 
		 * When dot product is zero, both signs are chosen to calculate multiply vertices.
		 * 
		 * eg. C 1,2 = Q(t) - 0.5*[sgnR(R.N) + sgnS(S.N) +/- T]
		 * 
		 * Repeat this for all 6 cube normals. If the case none return 2 dot products = 0, repeat with 
		 * all 12 edges, if the case none return  1 dot product, must be vertex collision.
		 * 
		 * Now include checks for finite plane.
		 * 
		 * For single vertex use standard Beck's algorithm that if Vertex.NL = 0 then true 
		 * where NL is each normal to the plane's multiply edges
		 * 
		 * If the Vertex !in Plane		=>			Plane < Cube
		 * Return Plane Vertex Closest to Cube Instead (For Face and Edge)
		 * 
		 * 
		 * For edge collision, both vertices must satisfy condition.
		 * 
		 * Disadvantage of method: If Face Size > Plane					(1)
		 * 							  Plane partially covers box		(2)
		 * 
		 * (2) can only occur for edge collision and face collision:
		 *
		 * For Edge collision, Check (normal to edge).(Plane Vertices) != 0.
		 * If this is false and the plane vertices that satisfy this condition are not = C 1,2, then include
		 * in collision points with box vertices involved in collision. Else ignore.
		 * 
		 * (2) impossible for face collision as face collision would not be detected.
		 * 
		 * (1) There are two solution to this:			Ignore
		 * 												Check for condition with C 1,2,3,4 calculation
		 * 
		 * Assuming you want to use an accurate physics engine (which is what this is!! kinda with some optimisation
		 * sacrifice) Then A completely new test condition must be raised. User can select where they want to 
		 * sacrifice optimisation for accuracy here??
		 * 
		 * 
		 *  Approximate Cube's CoM as sphere of radius r and linear vector Q(t). 
		 *  Perform Sphere to infinite plane L' collision test and return point of collision K. 
		 *  (if got here, means that collision with plane is possible!)
		 *  Find face normal, B, closest to collision point K along vector path.
		 *  For time t where K occurs, if B.(all plane vertices) is 0 & C1 2 3 4 dot product condition is false,
		 *  then face and face collision occurs and the plane is smaller than the face involved in collision.
		 *  Approximate plane as corner vertices and return all vertices.
		 *  
		 *  (1) & (2) hybrid scenario
		 *  
		 *  If Face collision detection is true for infinite plane but fails finite plane collision test, then finite plane
		 *  collision must be partially covering face. 
		 *  Find normal to each edge of plane and cube. Perform Beck's clipping using cube face edge normals to plane points.
		 *  Any points satisfy the condition was involved in collision.
		 *  Perform Beck's Clipping using plane's edge normals to cube's face vertices. Those that return true were involved
		 *  in collision. 
		 *  Return all points satisfying their respective conditions.
		 * 
		 */
		
		Maths T = new Maths();
		
		//Find normal to plane L
		//normal should be pointing towards initial position of box
		double[] Nl = Maths.Unit_Vector(Maths.Normal_to_Plane(Plane0));
		double[] CoM0 = T.CoM_Box_3D(Box0);
		double[] CoM1 = Maths.Vector_Add(CoM0,Box_Velocity);

		
		//normal should be pointing towards initial position of box
		//Move reference point to D
		//
		// D = normal*[x,y,z]
		//
		// CoM - N.D = Pact
		// N.Pact >= 0
		
		double[][] Plane1 = new double[Plane0.length][Plane0[0].length];
		int i;
		for(i = 0; i < Plane0.length; i++)
		{
			Plane1[i] = Maths.Vector_Add(Plane0[i],Plane_Velocity);
		}
		
		double D_Plane0 = ((Plane0[0][0]*Nl[0]) + (Plane0[0][1]*Nl[1]) + (Plane0[0][2]*Nl[2]));
		double D_Plane1 = ((Plane1[0][0]*Nl[0]) + (Plane1[0][1]*Nl[1]) + (Plane1[0][2]*Nl[2]));
		
		double[] Pact = Maths.Vector_Subtract(CoM0,Maths.Vector_Multiply(D_Plane0,Nl));
		double[] Nl_dot_CoM = Maths.Dot_Product(Pact,Nl);
		
		if( ((Nl_dot_CoM[0] < 0) || (Nl_dot_CoM[1] < 0)) || (Nl_dot_CoM[2] < 0))
			Nl = Maths.Vector_Multiply(-1,Nl);
		
		//Check CoM passes through L':		
		//
		//NL'.Q(0)/|NL'.Q(0)| != NL'.Q(1)/|NL'.Q(1)|
		//
		
		//Move reference point of Q to World Origin
		double[] Q_0 = Maths.Vector_Subtract(CoM0,Maths.Vector_Multiply(D_Plane0,Nl));
		double[] Q_1 = Maths.Vector_Subtract(CoM1,Maths.Vector_Multiply(D_Plane1,Nl));
		
		//Check If Box passes through plane
		//
		//Assume Linear Motion for now => Nl_0 = Nl_1
		//Find Whether Box intersects iwth plane. 
		//
		// If (Box CoM0.Plane0 N)/|(Box CoM0.Plane0 N)| == (Box CoM1.Plane1 N)/|(Box CoM1.Plane1 N)| 
		//
		// => Collision
		
		Boolean Box_Pass_Through_Plane = true;
		double[] test0 = new double[3]; 
		double[] test1 = new double[3];
		double[] Q_0_dot_Nl = Maths.Dot_Product(Q_0,Nl);
		double[] Q_1_dot_Nl = Maths.Dot_Product(Q_1,Nl);

		for(i = 0; i < 3; i++)
		{
			test0[i] = Maths.NaN_Checker(Q_0_dot_Nl[i]/Maths.Magnitude_Scalar(Q_0_dot_Nl[i]));
			test1[i] = Maths.NaN_Checker(Q_1_dot_Nl[i]/Maths.Magnitude_Scalar(Q_1_dot_Nl[i]));
		}
		
		if( ((test0[0] == test1[0]) && (test0[1] == test1[1])) && (test0[2] == test1[2]) )
			Box_Pass_Through_Plane = false;
		
		if(Box_Pass_Through_Plane == false)
		{
			System.out.println("Box Does Not Pass Through Plane");
			//Print_Vector("test0: ", test0);
			//Print_Vector("test1: ", test1);
			return null;
		}
		
		//Calculate circumscribed sphere radius, r 
		//r = (sqrt(3)/2) * (|R.N| + |S.N| + |T.N|)
		//N = normal to edge pointing away from cube
		double[] RST = {Maths.Magnitude_Scalar(CoM0[0]) - Maths.Magnitude_Scalar(Box0[0][0][0]), Maths.Magnitude_Scalar(CoM0[1]) - Maths.Magnitude_Scalar(Box0[0][0][1]),Maths.Magnitude_Scalar(CoM0[2]) - Maths.Magnitude_Scalar(Box0[0][0][2])};
		//Print_Vector("RST", RST);
		
		//Find all 6 plane normals
		double[][] Normal_to_Face = new double[6][3];
		for(i =0; i < 6; i ++)
		{
			Normal_to_Face[i] = Maths.Unit_Vector(Maths.Normal_to_Plane(Box0[i]));
			Maths.Print_Vector( ("N" + i), Normal_to_Face[i]);	
		}
		
		double circumscribed_radius = Maths.Magnitude_3D(RST);
		//System.out.println("r: " + circumscribed_radius);
					
		/*Find time L' intersects with CoM - Time of collision if radius taken as circumscribed
		 * 
		For moving planes									=>		D(t) = D0 + Vdt
		Here:														Vd = D1 - D0
																	N.Q(t) + D(t) - r = 0
																	N.Q0 + (Vq.t).N + D0 + Vd.t - r = 0
																	t(Vq.N + Vd) = r - D0 - N.P0
		Rearrange to solve for t:									t = | (r - |D0 - N.Q0| ) / ( |Vq.N + Vd| ) |
		*/
		//			Where Q(t) = CoM
		//				  D = Distance from plane to // plane through origin
		//
		// Note: HERE APPROX. OBJECT PROBLEM STARTS
		//
		double[] Vq = new double[3];
		Vq = Maths.Vector_Subtract(CoM1,CoM0);
		double Vd = D_Plane1 - D_Plane0;
		
		double numerator = circumscribed_radius - Maths.Magnitude_3D(Maths.Dot_Product(CoM0,Nl)) - D_Plane0;
		double numerator1 = (circumscribed_radius - Maths.Magnitude_3D(Maths.Dot_Product(CoM0,Nl)) + D_Plane0);
		
		if( numerator > numerator1)
			numerator = numerator1;
		
		double denominator = Maths.Magnitude_3D(Maths.Dot_Product(Vq,Nl)) + Vd;

		double t_collision = Maths.Magnitude_Scalar(numerator/denominator);

		//System.out.println("numerator: " + numerator);
		//System.out.println("denominator: " + denominator);
		System.out.println("t: " + t_collision);
		
		
		//Calculate Vertex Positions on collision
		
		/* If 0 <= t <= 1 then collision occurs. Else return null.
		 * 
		 * Zn = CoM(t) + [ +/- R +/- S +/- T]
		 * Where Zn - Vertex n position on collision
		 * 
		 * Check orientation of box on collision. Find signs of [R,S,T] such that NL'.Zn -> -oo
		 * If dot product > 0 => invert signs
		 */
		
		//Check condition if collision occurs within collision frame
		if( t_collision < 0 || t_collision > 1)
			return null;
		
		//Check Box orientation on collision
		//
		//RST.N must be negative. Check sign is positive. If true invert.
		// N is normal to cube's edge
		//
		//Perform Check
		//RST.Nl -> -oo
		//
		
		//Nl is pointing towards cube.
		//
		//Map RST to World (from reference of infinite plane @ D(t)
		//
		// D(t) = D0 + Vd.t
		
		double D_Collision = D_Plane0 + Vd*t_collision;
		
		//RST can only be positive real number
		//Map RST to world (from plane @t)
		//
		// RST_act = RST - N.D(t)
		//
		// NOTE: PLANE IS INFINITE: thus any point along collision face/edge can be calculated as RST of collision vertex
		//
		for(i = 0; i < 3; i++)
		{
			RST[i] = Maths.Magnitude_Scalar(RST[i]);
		} 			
		double[] RST_actual = Maths.Vector_Add(RST,Maths.Vector_Multiply(D_Collision,Nl));
		//Print_Vector("N.D(t)",Vector_Multiply(D_Collision,Nl));
		
		double[] Collision_RST = new double[3];
		for(i = 0; i < 3; i++)
		{
			if( (Nl[i]*RST[i]) > 0 )
				Collision_RST[i] = -RST[i];
			else
				Collision_RST[i] = RST[i];
		}
		
		//System.out.println("D(t): " + D_Collision);
		//Print_Vector("Nl",Nl);
		//Print_Vector("RST",RST);
		//Print_Vector("RST actual",RST_actual);
		//Print_Vector("Collision RST", Collision_RST);
		
		//Find CoM @ Collision
		//
		// Q(t) = CoM0 + Vcom.t
		// Vcom = CoM1 - CoM0
		double[] Vcom = Maths.Vector_Subtract(CoM1,CoM0);
		double[] CoM_Collision = Maths.Vector_Add(CoM0,Maths.Vector_Multiply(t_collision,Vcom));
		
		//Find normal to vertex
		//
		// (Zn - CoM_collision).N >= 0 if face normal is orthogonal to vertex
		//
		//Any vertex at a time can have up to a max of 3 normals
		//Vertex Normal must be dynamnically allocated
		//
		//Zn = CoM(t) + [ +/- R +/- S +/- T]
		
		double[] Zn = Maths.Vector_Add(CoM_Collision,Collision_RST);
		// Mapped_Zn_to_World = Collision_RST
		double[] Vertex_Normal = {0,0,0};
		double[] test;		
				
		//Nl is pointing towards CoM. Thus want opposite Nl.
		
		//Check Normals orthogonal to collision vertex is correct
		for(i = 0; i < 6; i++)
		{
			test = Maths.Dot_Product(Normal_to_Face[i],Collision_RST);
			
			if( ((test[0] >= 0) && (test[1] >= 0)) && (test[2] >= 0) )
			{
				//From Vertex normals find normal pointing towards Plane, N		
				//
				// Distance to plane = Nl.D(t)_Plane
				// As dot product vector and point is distance magnitude of point in direction vector.
				// For thus to work, point must be magnitude
				//
				// Where D(t) = |D(t)|
				//
				double[] Nl_dot_Dt = Maths.Vector_Multiply(Maths.Magnitude_Scalar(D_Collision),Maths.Vector_Multiply(-1,Nl));
				
				//Now Nl is normal of plane pointing towards box involved in collision
				//
				//Thus Nl.D(t) is the mapped collision point with infinite plane L (w. ref to infinite plane L'/CoM @origin)
				test = Maths.Dot_Product(Normal_to_Face[i], Nl_dot_Dt);
				
				//Normal with largest magnitude must be in direction towards C => collision vertex normal
				//As any vertex is connected by 3 faces => maximum of 3 calls
				if( Maths.Magnitude_3D(test) > Maths.Magnitude_3D(Vertex_Normal))
					Vertex_Normal = Normal_to_Face[i];
				
				//Print_Vector("Nl.D(t)",Nl_dot_Dt);
				//Print_Vector("N",Normal_to_Face[i]);
			}
		}	
		
					
		 /* Once [R,S,T] signs found, do dot product with normal to cube edge to find cube's orientation.
		 * 
		 * No dot product = 0 			=>			Single Vertex Collision
		 * 1 dot product = 0			=>			Edge Collision (Two Vertex)
		 * 2 Dot product = 0			=>			Face Collision ( 4 Vertex)
		 * 3 Dot Product = 0			=>			Plane is cube and aligned with Cube
		 * 
		 * Dot Products are |N.R|, |N.S| & |N.T|
		 * Where N is normal to cube's edge's face pointing away from CoM
		 */
		
		// Count number of vertices involved in collision.
		//
		// 2^No = Number of collision vertices
		// N.RST = dot product condition
		//
		// Where N = Collision Vertex Normal
		// Where No = Number of dot products = 0
		// 
		
		test = Maths.Dot_Product(RST,Vertex_Normal);
		int No_Dot_Products_0 = 0;
		
		//No for loop as doing long way avoid i increments computation
		if(test[0] == 0)
			No_Dot_Products_0++;
		if(test[1] == 0)
			No_Dot_Products_0++;
		if(test[2] == 0)
			No_Dot_Products_0++;
		
		//int Number_of_Collision_Vertices = (int)Power(2,No_Dot_Products_0);
		
		System.out.println("No_Dot_Products_0: " + No_Dot_Products_0);
		
		/* C = Q(t) - (sqrt(3)/2)*[sgnR(R.N) + sgnS(S.N) + sgnT(T.N)]
		 * 
		 * When dot product is zero, both signs are chosen to calculate multiply vertices.
		 * 
		 * eg. C 1,2 = Q(t) - (sqrt(3)/2)*[sgnR(R.N) + sgnS(S.N) +/- T]
		 */
		
		//Use number of vertices to decide collision array
		//
		// sgn = RST(t)/|RST(t)|
		//

		double[] sgn = new double[3];
		int[] element_dot_product_equal_zero = new int[No_Dot_Products_0]; 
		double[] sgnRSTxRST_dot_N = new double[3];
		
		//Calculate N.RST
		double[] coefficient = Maths.Dot_Product(Collision_RST,Vertex_Normal);
		int j = 0;		
		
		for(i = 0; i < 3; i++)
		{
			//Find sign of Collision RST			-		sign = RST/|RST|
			sgn[i] = Collision_RST[i]/Maths.Magnitude_Scalar(Collision_RST[i]); 
					
			//Is N.RST = 0
			if(coefficient[i] == 0)
			{
				//Then No N.RST. Value RST of C equation is RST
				sgnRSTxRST_dot_N[i] = RST[i];
				
				//If N.RST = 0, store element where it occurred.
				element_dot_product_equal_zero[j] = i;
				System.out.println("e" + j + ":" + element_dot_product_equal_zero[j]);
				j++;
			}
			// N.RST != 0
			else
			{
				//Find sign(RST).(N.RST)
				System.out.println(i + " is not dot product 0");
				sgnRSTxRST_dot_N[i] = sgn[i]*coefficient[i];
			}
			
		}
		

		double[][] Collision_Points_to_Infinite_Plane;
		
		//Calculate Infinite Plane Collision Points
		switch(No_Dot_Products_0)
		{
			case 8:
				Collision_Points_to_Infinite_Plane = new double[8][3];
				double[][] Vn = new double[8][3];
				for(int n = 0; n < 7; n++)
				{
					Vn[n][0] = CoM_Collision[0] + RST[0]*Maths.Power((-1),(n&1));
					Vn[n][1] = CoM_Collision[1] + RST[1]*Maths.Power((-1),(n&2));
					Vn[n][2] = CoM_Collision[2] + RST[2]*Maths.Power((-1),(n&4));
				}
				Collision_Points_to_Infinite_Plane = Vn;
				break;
			default:
				//C.Point = new double[(int) Power(2,No_Dot_Products_0)][3];
				Collision_Points_to_Infinite_Plane = Find_Edge_or_Face_Collision_Vertices(element_dot_product_equal_zero, CoM_Collision, sgnRSTxRST_dot_N);
				break;
		}
		
		//Print Collision Details for Debug
		if(Collision_Points_to_Infinite_Plane != null)
		{
			for(i = 0; i < Collision_Points_to_Infinite_Plane.length; i++)
			{
				Maths.Print_Vector( ("C" + i), Collision_Points_to_Infinite_Plane[i]);
			}
		}
		
		/* Now include checks for finite plane.
		 * 
		 * For single vertex use standard Beck's algorithm that if Vertex.NL = 0 then true 
		 * where NL is each normal to the plane's multiply edges
		 */
		
		/* Vertex - 2 Conditions:
		 * 
		 * 1 - Not Touching, Infinite Plane collision not finite 		=>		No Corners
		 * 2 - Touching													=> 		Cube Corner 
		 * 
		 * Edges & Faces - 4 Conditions: 
		 * 1 - Touching fully, Plane > Cube								=>		Cube Corners
		 * 2-				   Plane < Cube 							=> 		Plane Corners
		 * 2 - Touching Partially, Plane offset Cube					=> 		Plane & Cube Corners
		 * 4 - Not Touching, Touching Infinite Plane but not finite.	=>		No Corners
		 * 
		 * Vertices
		 */
		
		//Calculate Position of Plane @t
		//
		// P(t) = P0 + (P1-P0)*t
		
		double[][] Plane_at_Collision = new double[Plane0.length][3];		
		for(i  = 0; i < Plane0.length; i ++)
		{
			Plane_at_Collision[i] = Maths.Vector_Add(Plane0[i],Maths.Vector_Multiply(t_collision,Maths.Vector_Subtract(Plane1[i],Plane0[i])));
			Maths.Print_Vector(("Plane P" + i + " @t"),Plane_at_Collision[i]);
		}
		
		//Calculate Normal to each Plane @t edges
		//
		// Nl.En		where		En = Edge n of plane
		double[][] Plane_Edge_Normals = new double[Plane0.length][3];
		Plane_Edge_Normals = Maths.Calculate_Normal_to_Plane_Edges(Plane_at_Collision,Nl);
		
		
		//Handle Collision Conditions with finite plane
		double[] Vertex_ref_Plane_t;
		
		//For circumscribed approximation, Plane will be offset by circumscribed radius - RST
		//
		//To compensate for this;
		//
		// offset = (Dplane - Dcube).Nl
		// Where Nl is plane normal in direction of cube CoM
		//
		
		//Find Normal to Face involved in collision
		double[] Ncubeface = Maths.Unit_Vector(Maths.Normal_to_Plane(Collision_Points_to_Infinite_Plane));
		
		//Find D of Collision Face @t
		double D_cube_face_edge_on_collision = Maths.Magnitude_Scalar((Ncubeface[0]*Collision_Points_to_Infinite_Plane[0][0]) + (Ncubeface[1]*Collision_Points_to_Infinite_Plane[0][1]) +  (Ncubeface[2]*Collision_Points_to_Infinite_Plane[0][2]));
		D_Collision = Maths.Magnitude_Scalar(D_Collision);
		double[] Collision_Offset = Maths.Vector_Multiply((D_Collision - D_cube_face_edge_on_collision),Nl);
		
		System.out.println("Dcube: " + D_cube_face_edge_on_collision );
		Maths.Print_Vector("Collision Offset",Collision_Offset);
		
		//Offset Plane @t
		//
		//Note: Approximation happens here. Plane's origin vector path has been offset to cube due to overapproximation
		//
		for(i = 0; i < Plane0.length; i++)
		{
			Plane_at_Collision[i] = Maths.Vector_Add(Plane_at_Collision[i],Collision_Offset);
			Maths.Print_Vector(("offset P" + i),Plane_at_Collision[i]);
		}
		
		for(i = 0; i < Collision_Points_to_Infinite_Plane.length; i++)
		{
			System.out.println("----" + "Collision Pt(" + i + ")" + "-------");
			Maths.Print_Vector("Actual Collision Point",Collision_Points_to_Infinite_Plane[i]);
			
			//Check Conditions above 1 or 2
			//Map to world frame w. ref. plane
			//
			// Pact = C(t) - Nl.Dedge
			//
			// Due to offset, Dcube = Dplane
			//Apply Beck's Clipping Algorithm
			//
			// Nedge.C(t) <= 0
			//
			Boolean Collision = true;

			for(j = 0; j < Plane0.length; j++)
			{
				double Dedge = (Plane_at_Collision[j][0]*Plane_Edge_Normals[j][0]) + (Plane_at_Collision[j][1]*Plane_Edge_Normals[j][1]) + (Plane_at_Collision[j][2]*Plane_Edge_Normals[j][2]);
				double[] Dedge_dot_Nl = Maths.Vector_Multiply(Dedge,Nl);
				Pact = Maths.Vector_Subtract(Collision_Points_to_Infinite_Plane[i],Pact);					
				//System.out.println("Dedge: " + Dedge);
				Maths.Print_Vector("Pact",Pact);
				
				test = Maths.Dot_Product(Pact,Plane_Edge_Normals[j]);
				
				if( (test[0] < 0 || test[1] < 0) || test[2] < 0)
				{
					Collision = false;
					Maths.Print_Vector("Plane Nornal (No C)", Plane_Edge_Normals[j]);
					Maths.Print_Vector("C (No C)",Collision_Points_to_Infinite_Plane[i]);
				}
					
			}
			
			switch(Collision_Points_to_Infinite_Plane.length)
			{
			//Vertex Case - Either Collided Or Not.
			//Beck's Clipping with reference to Plane Edge involved in calculation
			//
			case 1:
				
				break;
				
			//Edge/Face Collision Conditions
		/*	case 2: case 4:
				
				for(j = 0; j < Plane0.length; j++)
				{
					//D_Edge = (Plane_Edge_Normals[j][0]*Plane_at_Collision[j][0]) + (Plane_Edge_Normals[j][1]*Plane_at_Collision[j][1]) + (Plane_Edge_Normals[j][2]*Plane_at_Collision[j][2]);
					//Vertex_ref_Plane_t = Vector_Subtract(Collision_Points_to_Infinite_Plane[i],Vector_Multiply(D_Edge,Nl));
					
					//Compensate for radius offset with plane
					//
					// Cube Vertex @t - Plane Collision Point
					//double[] approximation_offset_at_collision = Vector_Subtract()
					
					System.out.println("D @t:" + D_Collision);
							
					Vertex_ref_Plane_t = Vector_Subtract(Collision_Points_to_Infinite_Plane[i],Plane_at_Collision[j]);
					
					//System.out.println("D_Edge: " + D_Edge);
					Print_Vector("World Collision Pt",Vertex_ref_Plane_t);
					test = Dot_Product(Plane_Edge_Normals[j],Vertex_ref_Plane_t);
					//Print_Vector("test:",test);
					
					//N.Vertex >= 0 for collision
					if( ((test[0] < 0) || (test[1] < 0)) || (test[2] < 0) )
					{
						System.out.println("No Finite Plane Collision");
						//Print_Vector("Failed Case",test);
					}
					
				}

				//System.out.println("Vertex " + i);
				//Print_Vector("Nl.D",Vector_Multiply(D_Collision,Nl));
				//Print_Vector("Vref",Vertex_ref_Plane_t);
				*/
			default:
					break;
			}
		}
		
		return C;
	}
	
	//Returns vertices where dot product = 0 is +/- RST
	//
	//Inputs: Element of N.RST = 0, Q(t), sgnRST(N.RST_Collision),
	
	static double[][] Find_Edge_or_Face_Collision_Vertices(int[] Element_equal_0, double[] CoM, double[] RST)
	{
		
		Maths T = new Maths();
		
		int Number_of_Collision_Vertices = (int)Maths.Power(2,Element_equal_0.length);
		//System.out.println("Number of Collision Vertices: " + Number_of_Collision_Vertices);
		//System.out.println("Element_equal_0.length: " + Element_equal_0.length);
		double[][] Point = new double[Number_of_Collision_Vertices][3];
		int j = 0;
		
		switch(Number_of_Collision_Vertices)
		{
			//Vertex Collision
			case 1:
				Point[0] = Maths.Vector_Add(CoM,RST);
				break;
			//Edge Collision
			case 2:
				Point[j] = Maths.Vector_Add(CoM,RST);
				RST[Element_equal_0[0]] = -RST[Element_equal_0[0]];
				Point[++j] = Maths.Vector_Add(CoM,RST);
				break;
			//Face Collision
			case 4:
				//Find element that is constant
				//
				// 3 - |Sum of dot product = 0 elements| = constant element
				int Constant_Element = (int) (3 - Maths.Magnitude_Scalar(Element_equal_0[0] + Element_equal_0[1]));
				System.out.println("constant element: " + Constant_Element);
				
				// Use bitwise operator over 2 variables to calculate collision points
				// -1^(n & 1) && -1^(n&2)
				for(int n = 0; n < 4; n++)
				{
					Point[n][Element_equal_0[0]] = CoM[Element_equal_0[0]] + (RST[Element_equal_0[0]]*Maths.Power(-1,(n&1)));
					Point[n][Element_equal_0[1]] = CoM[Element_equal_0[1]] + (RST[Element_equal_0[1]]*Maths.Power(-1,(n&2)));
					Point[n][Constant_Element] = CoM[Constant_Element] + RST[Constant_Element];
				}
				break;
			default:
				break;
		}
		
		return Point;
	}
	
	//For moving box calculates point of collision.
	//
	//Inputs: Sphere start, Sphere end, box start, box end, radius of sphere
	//
	//Return: Class containing Collision Data
	
	//REGRESSION TEST TO SEE IF THIS ACTUALY WORKS!!!! - MAY BE MIRACLE FASTER
	// Becasue if collided with any plane, must be a box collision!
	
/*	static Collision_Information Sphere_to_Moving_Box(double[] P0, double[] P1, double[][][] C0, double[][][] C1, double r)
	{
		//Find time at which collision occurs with box surface
		double t = Find_Time_Moving_Box_Sphere_Collides(P0,P1,C0,C1,0);
		
		//Calculate point, P, at this time of collision
		double[] Vp = Vector_Subtract(P1,P0);
		double[] P = Vector_Add(P0, Vector_Multiply(t,Vp));
		
		//Store Collision point and time in Collision Class
		Collision_Information C = new Collision_Information();
		C.Point = P;
		C.time = t;
		
		return C;
	} */
	
	static Collision_Information Linear_Sphere_to_Moving_Box(double[] P0, double[] P_Velocity, double[][][] C0, double[] C_Velocity, double r)
	{
		
		Maths T = new Maths();
		
		//Find time at which collision occurs with box surface
		double t = Find_Time_Moving_Linear_Box_Sphere_Collides(P0,P_Velocity,C0,C_Velocity,0);
		
		//Check collision is within frame
		// 0 <= t <= 1
		if( t < 0 || t > 1)
			return null;
		
		//Calculate Box Linear Velocity
		//Assume Translation Matrix constant for all points
		
		double[] Vc = C_Velocity;
		
		//Translate Box to position when collision occurred.
		double[][][] C_at_Collision = new double[6][4][3];
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 4; j++)
			{
					C_at_Collision[i][j] = Maths.Vector_Add(C0[i][j],Maths.Vector_Multiply(t,Vc));
			}
		}
		
		//Find Collision details by doing Sphere to stationary box detection
		return Sphere_to_Stationary_Box(C_at_Collision, P0, P_Velocity, r);
	} 
	
	//For moving cube, maps point to world frame along box's velocity vector.
	//Calculates time at which Beck's clipping algorithm is true for a face along vector path
	//
	// C(t) = cube vertex = C0 + Vct
	// D = N.C(t)x + N.C(t)y + C(t).z
	// Pact = P(t) - N.D
	// Pact.N = 0 					-		when surface collision occurs as first pt. of collision.
	// Condition:	0 <= t <= 1
	// t = (N3.C0 - N.P0)/N.(Vp-VcN2)
	// L': N.P + D = r;
	//
	//Inputs: Sphere Start, Sphere End, Box Start, Box End, Radius
	
	static double Find_Time_Moving_Linear_Box_Sphere_Collides(double[] P0, double[] P_Velocity, double[][][] C0, double[] C_Velocity, double r)
	{
		
		Maths T = new Maths();
		
		//Calculate Vp (Sphere vector)
		double[] Vp = P_Velocity;
		
		//Calculate Translation Matrix for Box
		//Assume Translation Matrix is identical for all Box points
		double[] Vc = C_Velocity;
		
		//Ensure Normal is pointing away from CoM.
		double[] CoM = T.CoM_Box_3D(C0);
		double[][] Normal_of_Box = new double[6][3];
		double t = -0.000001;
		double lowest_t = t;
		
		for(int i = 0; i < 6; i++)
		{
			//Find Normal to Box (initial Box is fine as only deals with linear motion for time being)
			Normal_of_Box[i] = Maths.Unit_Vector(Maths.Normal_to_Plane(C0[i]));
			
			//Map CoM to Plane's local frame (D)
			//D = Normal*[x,y,z]
			//Pact = CoM - N.P
			double D = (Normal_of_Box[i][0]*C0[i][0][0]) + (Normal_of_Box[i][1]*C0[i][0][1]) + (Normal_of_Box[i][2]*C0[i][0][2]);
			double[] Pact = Maths.Vector_Subtract(CoM, Maths.Vector_Multiply(D,Normal_of_Box[i]));
		
			//If Pact.N < 0			=>			Normal pointing away from Normal

			double[] test_vector = Maths.Dot_Product(Pact,Normal_of_Box[i]);
			
			if( ((test_vector[0] > 0) || (test_vector[1] > 0)) || (test_vector[2] > 0) )
				Normal_of_Box[i] = Maths.Vector_Multiply(-1,Normal_of_Box[i]);
			
			
			//Find lowest value of t involved in collision
			//t = (N3.C0 - N.P0)/N.(Vp-VcN2)
			//As normal is normalized, any odd powers can be treated as N.
			double[] N3_dot_C0 = Maths.Dot_Product(Normal_of_Box[i],C0[i][0]);
			double[] N_dot_P0 = Maths.Dot_Product(Normal_of_Box[i],P0);
			double[] Vc_dot_N2 = Maths.Dot_Product(Vc,Maths.Dot_Product(Normal_of_Box[i],Normal_of_Box[i]));
			
			lowest_t = t;
			t = (Maths.Magnitude_3D(Maths.Vector_Subtract(N3_dot_C0,N_dot_P0)) - r) / Maths.Magnitude_3D(Maths.Dot_Product(Normal_of_Box[i],Maths.Vector_Subtract(Vp,Vc_dot_N2)));
		
			//Note:
			//If t returns infinity, likely that Box/Point has no components in a plane
			
			//Find point of first collision: 
			//when t -> 0  but t >= 1
			if(t != Double.POSITIVE_INFINITY && t != Double.NEGATIVE_INFINITY)
			{
				if( (lowest_t > t) && (t >= 0) )
					lowest_t = t;
			}
			
			/*System.out.println("----");
			Print_Vector( ("n" + i), Normal_of_Box[i]);
			Print_Vector("N3.C0: ", N3_dot_C0);
			Print_Vector("N.P0: ", N_dot_P0);
			Print_Vector("Vc.N2: ", Vc_dot_N2);
			System.out.println("t: " + t);*/
		}
		
		System.out.println("----");
		System.out.println("t_lowest: " + lowest_t);
		return lowest_t;
	}
	
	// Calculates Collision point between Sphere and Box
	//
	//Inputs:	Box Vertex Information 3D Double	[face][vertex][x,y,z component]	  -	Assumes Cube is already sorted in this format
	//			Sphere Initial Position 3D double (x,y,z)
	//			Sphere Final Position 3D Double (x,y,z)
	//			Sphere Radius
	//
	//Returns:	Point of Collision		-	double[3] = (-1,0,0) if no collision
	//
	//Return:	If Collision		-			Point of collision (x,y,z) & time of collision (t) as double[4] respectively
	
	static Collision_Information Sphere_to_Stationary_Box(double[][][] Box, double[] Sphere_Initial, double[] Sphere_Velocity, double r)
	{
		System.out.println("----");
		
		Collision_Information C_Data = new Collision_Information(1);
		double[][] cube_normal = new double[6][3];
		double[] n = new double[3];
		double[] P1 = new double[3];
		double[] P0 = new double[3];
		
		//If Booleans change value			=>			Collision has occurred
		Boolean collision_occurred = true;
		Boolean point_passed_through_box = false;
		
		System.out.println("----");
		
		//Beck Clipping to calculate whether point at t=0 & t=1 has collided with Box
		//
		//First calculate normal of each cube's face
		//
		//Two lines involved in normal calculation must share a vertex.
		//
		//Connect the shared vertex with any two other vertex. 
		//Two lines, A and B, parallel to any of the face's edges should be generated.
		//These conditions must stand if the normal calculated is perpendicular to the face and pointing away from the cube.
		//
		//	A.B = 0 			=>				Two lines will generate surface normal perpendicular to face (line a || b is across face not aligned)
		//  (A x B) = N
		//	N.P < 0				=>				Normal pointing away from Box's CoM
		//
		//  P = any point on the opposite face in the direction -N [In this algorithm, a vertex is used as P]
		//
		// Repeat for all 6 faces to get 6 normals.
		//
		// Preliminary Step - Use Beck Clipping algorithm to determine if P1 is inside of Cube. Assume no collision at P0.
		//
		// for(x=1:6) { (N of face x).P1 <= r) }		-		All conditions should stand if collision occurs
		//
		// r = sphere's radius
		//
		// if <= r, dot product -> 0 is plane of collision
		//
		// Any condition fails at t = 0 or t = 1			=>			Do Binary Search to see if sphere passed completely through box or doesn't touch box at all
		//
		// Do binary search for 3 or 5 iterations using above conditions for each new value of t between 0 and 1.
		// If above conditions doesn't stand for all values of t		=>		No Collisions.
		// Else if conditions stands at one t							=> Use normal @this t whose N.P(t) is closest to 0 in plane calculation [only negative -> 0 values are valid]
		// Else if several t's stand 									=> Use normal whose N.P(t) is closest to 0 (at any t) in plane calculation [only negative -> 0 values are valid]
		//
		// The normal whose N.P(t) is closest to 0 is the first surface involved in collision. [only negative -> 0 values are valid]
		// For this face, use it's vectors A and B to generate plane L = {N, D}.
		// This is collision plane.
		// Do 2D plane to sphere calculation on this plane:
		//
		// L = {N,D}
		// L' = {N, (D - r)}
		// L' => N.P(t) + D - r = 0
		//		 N.P0 + Vt.N + D - r
		//		 t = | (r - D - N.P0) / N.V |
		//	C = P(t) + r.N
		//
		// C is point of intersection
		// t is time of intersection along sphere's V
		
		// Calculate Cube's CoM
		//
		// CoM = [ (Ax - Bx), (Ay - Cy), (Az - Dz) ] 
		//
		// A,B,C,D represent vertices of the box where same plane components are not equal (eg. x of A != x of B)
		//
		// Maximum of 2 faces are required to get change in x, y and z. This holds true for all signs of the two vertex signs.
		//
		
		Maths T = new Maths();
		
		double[] CoM = new double[3];
		CoM = T.CoM_Box_3D(Box);
		
		Maths.Print_Vector("CoM", CoM);
		
		double[] a = new double[3];									// Define all lines making plane, normals and opposite face's edge and initialise all 1,1,1 to purposesly not obey below conditions
		double[] b = new double[3];
		double[] Sphere_Final;
		
		// Find each normal to cube pointing away from cube's centre of mass (CoM)
		for(int i = 0; i < 6; i++)										// Find normals pointing outwards for all 6 faces
		{
			int j = 1;

			a[0] = b[0] = a[1] = b[1] = a[2] = b[2] = 1;
			
			// Create two lines defining plane and check they are perpendicular - Max 3 loops
			n = Maths.Normal_to_Plane(Box[i]);
				
				//P.N <= 0 if N is pointing out of cube. P is cube's CoM. 
				//CoM is translated by face's vertex to keep same reference point with normal
				if(Maths.Magnitude_3D(Maths.Dot_Product(n,(T.Vector_Subtract(CoM,Box[i][0])))) > 0)					// If N.P > 0, N must be pointing towards cube's centre
					n = Maths.Vector_Multiply(-1,n);
				
				//Assign n to corresponding face.
				n = Maths.Unit_Vector(n);
				cube_normal[i][0] = n[0]; 
				cube_normal[i][1] = n[1]; 
				cube_normal[i][2] = n[2]; 
				
	//		System.out.println("i: " + i);
			Maths.Print_Vector(("n" + i), n);
	//		Print_Vector("cube_normal: ", cube_normal[i]);
				
				//Find whether point has intersected with cube - 
				//take point reference at ANY edge of the plane involved in calculation
			
				//Initialise new set of P(0) & P(1) with reference to box face vertex
				
			Sphere_Final = Maths.Vector_Add(Sphere_Initial, Sphere_Velocity);
			
				P1 = Maths.Vector_Subtract(Sphere_Final,Box[i][0]);
				P0 = Maths.Vector_Subtract(Sphere_Initial,Box[i][0]);				
				
				//Beck Clipping Condition: P(1).N <= r 			=>				Collision
				double[] test_vector_P1 = Maths.Dot_Product(n,P1);
				double[] test_vector_P0 = Maths.Dot_Product(n,P0);
				
				if( ((test_vector_P1[0] > r) || (test_vector_P1[1] > r)) || (test_vector_P1[2] > r))
				{
					//Print_Vector("N.P1", test_vector_P1);
					Maths.Print_Vector("No Collision normal", n);
					collision_occurred = false;
				}
				
				//Preliminary Check
				//If no collision occurred, check point did not pass through cube.
				//P(0).N/|P(0).N| == P(1).N/|P(1).N| 			=>			No Collision as stays on same side of bounding box
				
				//initialise variables
				double[] sign_P1 = new double[3];
				double[] sign_P0 = new double[3];
				test_vector_P1 = Maths.Dot_Product(n,test_vector_P1);
				test_vector_P0 = Maths.Dot_Product(n,test_vector_P0);
				
			//	Print_Vector("N.P1", test_vector_P1);
			//	Print_Vector("N.P0", test_vector_P0);
				
				//Calculate sign
				for(int k = 0; k < 3; k++)
				{
					sign_P1[k] = Maths.NaN_Checker(Maths.Sign_Magnitude_Scalar(test_vector_P1[k]));
					sign_P0[k] = Maths.NaN_Checker(Maths.Sign_Magnitude_Scalar(test_vector_P0[k]));
				}
				
				//stdout debug
			//	Print_Vector("sign_P0", sign_P0);
			//	Print_Vector("P0", P0);
			//	Print_Vector("sign_P1", sign_P1);
			//	Print_Vector("P1", P1);
				
				//Condition Check
				if( ((sign_P1[0] != sign_P0[0]) || (sign_P1[1] != sign_P0[1])) || ((sign_P1[2] != sign_P0[2])))
				{
					System.out.println("Point Passed Through Box");
					point_passed_through_box = true;
				}
		}

		//Preliminary Check Condition
		if(collision_occurred == false && point_passed_through_box == false)
			return null;
		
		System.out.println("----");
					
		//If got to here => collision has occurred
		//If collision occurred, for time t, find plane closest to collision point - N.P1 -> 0   - Binary Search
		//Binary Search iteration accuracy set by user (no sig fig)
		//
		//Collision Frame always between 0 <= t <= 1
		//
		//P(t) = P(0) + Vp.t			where			Vp = P(1) - P(0)
		double[] Vp = Sphere_Velocity;
		int face_involved_in_surface_collision = binary_search_box(0, 1, cube_normal, Box, Sphere_Initial, Vp);
		int i = face_involved_in_surface_collision;
		System.out.println("index of face closest to collision: " + face_involved_in_surface_collision);
		
		Plane_Handler H = new Plane_Handler();
		
		double[] Plane_Velocity = {0,0,0};
		
		//For plane closest to point at collision, do sphere -> 2D plane collision check
		C_Data = Plane_Handler.Linear_Sphere_to_3D_Finite_Plane(Sphere_Initial, Sphere_Velocity, Box[i], Plane_Velocity, r);
		
		return C_Data;
	}
	
	
	//Finds the face a sphere intersects with before interpenetrating a bounding box.
	// Percision (sig fig) set by user.
	//
	//Returns index face corresponds to
	//
	//If returns -1, (Should Never Due to Our Box Collision Conditions!!!!!), cannot find collision plane
	
	static int binary_search_box(double imin, double imax, double[][] normal, double[][][] edge, double[] P0, double[] Vp)
	{
		
		Maths T = new Maths();
		
	  // continue searching while [imin,imax] is not empty
		double lowest_dot_product = 100;
		int face_number_closest_to_collision_point = -1;
		
		int binary_search_iteration = 0;
		double[] P_t = new double[3];
	    double D; double[] N_dot_D = new double[3];
		
		
	//Collision frame always occurs between 0 <= t <= 1
	  while (imax >= imin)
	    {
		  binary_search_iteration++;
		  
	      /* calculate the midpoint for roughly equal partition */
	      double t = ((imax-imin)/2) + imin;
	      
	      //Store previous dot product for tmax
	      double previous_dot_product = lowest_dot_product;
	      
	      //For P(t), Find lowest dot product N.P(t) for each face
	      for(int i = 0; i < 6; i++)
	      {
		      // determine which subarray to search
		      //P(t) = P(0) + Vp.t with reference to face
		      P_t = Maths.Vector_Add(P0,Maths.Vector_Multiply(t,Vp));
	    	  
		      //D = N*[x,y,z]			then				P(t) = P(t) - N.D			to get reference point
		      //
		      //N.D = distance between face and // plane to face passing through origin
		      //
		      //|P(t) - N.D| -> 0
		      D = (edge[i][0][0]*normal[i][0]) + (edge[i][0][1]*normal[i][1]) + (edge[i][0][2]*normal[i][2]);
		      N_dot_D = Maths.Vector_Multiply(D,normal[i]);
		      P_t = Maths.Vector_Subtract(P_t,N_dot_D);
		      
		      //For P(t), Find lowest dot product N.P(t) for each face
	    	  double dot_product_magnitude = Maths.Magnitude_3D(P_t);
	    	  if(dot_product_magnitude < previous_dot_product)
	    	  {
	    		  lowest_dot_product = dot_product_magnitude;
	    		  face_number_closest_to_collision_point = i;
	    	  }
	    	 
	      }
	      
	      //If P(t).N > P(t_previous).N, then collision occurs between tmid -> tmax.
	      //else collision occurs between tmin -> tmid
	      //
	      // Allow N.P(t) -> 0 precision to 5 sig figures
	      
	      //Round dot product result to 5 sig fig
	      BigDecimal bd = new BigDecimal(previous_dot_product);
	      bd = bd.round(new MathContext(5));
	      previous_dot_product = bd.doubleValue();
	      BigDecimal b = new BigDecimal(lowest_dot_product);
	      b = b.round(new MathContext(5));
	      lowest_dot_product = b.doubleValue();
	      		      
  //  	  System.out.println("*****");
   // 	  System.out.println("t: " + t);
   // 	  Print_Vector("P(t)",P_t);
//			System.out.println("lowest N.P(t): " + lowest_dot_product);
//	      System.out.println("previous N.P: " + previous_dot_product);
	      
	      if      ( lowest_dot_product <  previous_dot_product)
	      {
	        // change min index to search upper subarray
	        imin = t;
	      }
	      else if ( lowest_dot_product > previous_dot_product )
	      {
	        // change max index to search lower subarray
	        imax = t;
	      }
	      else
	      {
	    	  System.out.println("binary_search_iteration: " + binary_search_iteration);
	    	  // key found at index imid
		        return face_number_closest_to_collision_point;
	      }
	   
	    }
	  // key not found
//	  System.out.println("binary_search_iteration: " + binary_search_iteration);
	  return -1;
	}
	
	
}
