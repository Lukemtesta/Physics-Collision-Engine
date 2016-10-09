
//--------------------------------------------------Updated Sphere-to-Sphere 24.2.13------------------------------------------

public class Sphere_Handler {
		
			 /* Bounding Sphere: Y
			 * 
			 * Sphere to Sphere Y (Main)
			 * Sphere to Point Y
			 * Point to Point Y
			 * Sphere to 1D Plane Y
			 * Point to 1D Plane Y
			 * 1D plane to 1D plane Y
			 * Sphere to 2D Plane Y (Main)
			 * Point to 2D Plane Y
			 */
	
	
			//Calculates Collision point between Two Spheres - Stationary or Moving
			
			//Inputs:   Sphere 1 first and final position as array respectively (3 inputs: x,y,z) - P(t)
			//			Sphere 2 first and final position as array respectively (3 inputs: x,y,z) - Q(t)
			//			Sphere 1 radius
			//			Sphere 2 radius
			//
			//Return:	double[3] = (-1,0,0) if no collision
			//
			//Return:	If Collision		-			Point of collision (x,y,z) & time of collision (t) as double[4] respectively
						
			static Collision_Information Sphere_to_Sphere(double[] Sphere_1_Initial, double[] Sphere_1_Velocity, double[] Sphere_2_Initial, double[] Sphere_2_Velocity, double R1, double R2)
			{
				
				Maths T = new Maths();
				
				//P(t) is Linear Vector Describing Sphere 1's movement from Sphere 1 centre of mass;
				//P(0) = Initial position
				//P(1) = Final Position
				//P(t) = P0 + Vp*t	where Vp = P1 - P0
				
				//Calculate Gradient Vp. Stored as double x, y, z;
				double[] Vp = Sphere_1_Velocity;
				
				//Q(t) is Linear Vector Describing Sphere 2's movement from Sphere 2 centre of mass
				//Q(0) = Initial position
				//Q(1) = Final Position
				//Q(t) = Q0 + Vq*t	where Vq = Q1 - Q0
					
				//Calculate Gradient Vq. Stored as double x, y, z;
				double[] Vq = Sphere_2_Velocity;
				
				// Collision Distance Condition, d = R1 + R2
				// ||R1 + R2||^2 = R_Condition_Squared
				
				double R_Condition = R1 + R2;
				double R_Condition_Squared = R_Condition*R_Condition;
				
				//D^2 = || P(t) - Q(t) ||^2
				//D^2 = A^2 + 2t*(A.B) + (t^2)*(B^2)
				//B = Vp - Vq
				//A = P0 - Q0
				
				//Calculate B
				double[] B = T.Vector_Subtract(Vp,Vq);
				
				//Calculate A
				double[] A = T.Vector_Subtract(Sphere_1_Initial,Sphere_2_Initial);
				
				//Solving for t: D^2 = A^2 + 2t*(A.B) + (t^2)*(B^2)
				//		
				//			=> t1 = [ -(A.B) - sqrt( (A.B)^2 - (B^2) * ( A^2 - d^2) ) ] / B^2		- Using Quadratic Formula
				//
				//t1 = towards motion collision tangent
				//t2 = away motion collision tangent
				//t1 = t2 single surface collision
				//t1 != t2 two collision points (interpenetration occurred)
				//
				// => Only t1 Interested (First Point of Collision)
				//
				//	0 < t1 < 1					=> Collision
				//  t1 < 0 || t1 > 1			=> No Collision
				//
				//	B^2 > 0						=> t1 <= t2
				//	B^2 = 0 					=> No Collision - Stationary or parallel motion
				//
				//	sqrt(input): input <= 0		=> No Collision
				
				
				//Derivative of D^2 computes time at which D^2 is Min - Preliminary Check
				//
				// d(D^2)/dt = 0 = 2(A.B) + 2t(B^2)
				//
				// solve t:		t = -(A.B) / (B^2)
				//
				// where t is time at which D^2 is minimum
				//
				// Dmin^2 <= |R1 + R2|^2		=>	Collision
				
				//Compute Required Variables for t_min_when_Distance_Between_Sphere_CoM_Squared_is_min	
				//
				//All variables are variation of A & B
				double[] B_squared = T.Dot_Product(B,B);	
				double[] A_dot_B = T.Dot_Product(A,B);
				double[] negative_A_dot_B = T.Vector_Multiply(-1,A_dot_B);
				
				//Find time at which Distance between sphere's CoM is minimum
				double tmin = T.Magnitude_3D(negative_A_dot_B)/T.Magnitude_3D(B_squared);
				
				// Calculate Dmin^2 from t @ Dmin
				//
				// Dmin = |P(t @ tmin) - Q(t @ tmin)|			<-------------
				//
				// Dmin_squared is used instead of Dmin because for Dmin_squared the below solution is valid
				//
				// Dmin^2 = |P(tmin) - Q(min)|^2
				//        = P0 + Vp.tmin - Q0 - Vq.tmin
				//        = (A + t.B)^2
				// f(t)   = A^2 + 2tmin.(A.B) + (t^2).(B^2)
				//
				// Then df(t)/dt = 0 for t when Dmin^2 is minimum
				//
				// df(t)/d(t) = 0 = 2.(A.B) + 2tmin.(B^2)
				//
				// Thus tmin = -(A.B)/B^2
				//
				// Alternative equation for Dmin in terms of A & B achieved by subbing tmin into the quadratic produces
				//
				// D = |P(t) - Q(t)| = | P0 - Q0 + Vp.t - Vq.t|
				//   = A + Bt
				// A = P0 - Q0
				// B = Vp - Vq
				//
				// D = A - (A.B_squared)/B_squared				<-------------
				//
				
				double[] P_at_tmin = T.Vector_Add(Sphere_1_Initial,T.Vector_Multiply(tmin,Vp));
				double[] Q_at_tmin = T.Vector_Add(Sphere_2_Initial,T.Vector_Multiply(tmin,Vq));
						
				double Dmin = T.Magnitude_3D(T.Vector_Subtract(P_at_tmin,Q_at_tmin));
				double D_squared = (R1 + R2)*(R1 + R2);
				
				System.out.println("D^2: " + D_squared);
				
				//Preliminary Check
				//
				// Check if Dmin results in collision along linear velocities
				if(Dmin > R_Condition)
				{
					System.out.println("No Collision,       Dmin^2 > |R1 + R2|^2:  " + Dmin + " > " + R_Condition_Squared + "                " + "tmin: " + tmin);
					return null;
				}
				else
				{
					System.out.println("Collision,   Dmin^2 <= |R1 + R2|^2:  " + Dmin + " <= " + R_Condition_Squared + "                " + "tmin: " + tmin);
				}
				
				//Compute Required Variables for t_collision	
				//
				//All variables are variation of A & B
				double[] A_squared = T.Dot_Product(A,A);
				double[] A_dot_B_squared = T.Dot_Product(A_dot_B,A_dot_B);
				
				//If it makes it this far. Collision occurs.
				
				//Calculate t at which surface collision occurs
				//
				// surface collision occurs when D = |rp + rq|
				//
				// Deriving t:
				//
				// D_squared = A_squared + 2t.A.B + B_squared.t_squared
				//
				// D_squared = (Rq + Rp)^2
				//
				// Sub D_squared into quadratic equation defining t to find collision time when surfaces collide (time @ D = Rq + Rp)
				//
				// (Rq + Rp)^2 = A_squared + 2t.A.B + B_squared.t_squared
				//
				// 0 = A_squared - (Rq + Rp)^2 + 2t.A.B + B_squared.t_squared
				//
				// Using Quadratic Formula at^2 + bt + c = 0
				//
				// a = B_squared
				// b = 2A.B
				// c = A_squared - (Rq + Rp)^2
				//
				// t @ D = ( -2A.B +/- Sqrt( 4[(A.B)^2] - 4(B^2)[ (A^2) - (Rq + Rp)^2 ] ) / B^2 
				//
				// t @ D = ( -2A.B +/- Sqrt( 4[(A.B)^2] - 4(B^2)[ (A^2) - (Rq + Rp)^2 ] ) / B^2
				//
				// Note (A.B)^2 = (A^2).(B^2)
				//
				// IMPORTANT: |A| = A & |B| = B else IT WILL NOT WORK!!!!
				//
				// If B^2 = 0	=>	At rest or moving parallel thus no collision
				//    B^2 > 0	=>  t1 < t2				=>			-b^2 - sqrt(coefficient)
				//
				// This B^2 condition works as B = Vp - Vq describing P & Q's motion
				
				double mag_A = T.Magnitude_3D(A);
				double mag_B = T.Magnitude_3D(B);
				
				double mag_A_squared_minus_D_squared = (mag_A*mag_A) - R_Condition_Squared;
				double mag_B_squared_dot_A_squared_minus_D_squared = (mag_B*mag_B) * mag_A_squared_minus_D_squared;
				double mag_A_dot_B_squared = (mag_A*mag_B)*(mag_A*mag_B);
				
				double sqrt_coefficient = mag_A_dot_B_squared - mag_B_squared_dot_A_squared_minus_D_squared;
				double t_collision = ( (mag_A*mag_B) - Math.sqrt(sqrt_coefficient) ) / (mag_B*mag_B);
				
				//Preliminary Check
				//
				// Sqrt must be real ( >= 0 )
				if(sqrt_coefficient < 0)
				{
					System.out.println("Collision Invalid, sqrt_coefficient < 0: " + sqrt_coefficient);
				}
				else
				{
					System.out.println("Collision True, sqrt_coefficient >= 0: " + "    " +  sqrt_coefficient + " >= 0");
				}
								
				//For Non-Constant Velocities, Time of collision must be 0<t<1 to have collided in animation frame.	
				//
				// Collision in animation frame
				//
				// 0 <= t <= 1
				//
				if( (t_collision < 0) || (t_collision > 1))
				{
					System.out.println("Collision Time Incorrent,   t < 0 || t > 1: " + "   " + t_collision + "< 0 || " + t_collision + "> 1");
					return null;
				}
				else
				{
					System.out.println("Collision Occurred,   0 <= t <= 1: " + "   0 <= " + t_collision + " <= 1");
				}
					
				//Once Collision Time, Find Q(t_collision) & P(t_collision)
				//Find Mid-Points of Sphere at collision time
				//
				// Q(t_collision) = Qc
				// P(t_collision) = Pc
				//
				// Collision point on vector connecting P(t) and Q(t) [taken as unit vector in direction towards P]
				//
				// L = P(t) - Q(t) / |P(t) - Q(t)| 		- where t = t_collision
				// 
				//
				// Where:	C = P(t_collision) + R1*L = Q(t_collision) + R2*L
				
				//Find P(Collision) & Q(Collision)
				//
				// P(t_collision) = P0 + Vp*t_collision
				//
				double[] P = T.Vector_Add(Sphere_1_Initial, T.Vector_Multiply(t_collision,Vp));
				double[] Q = T.Vector_Add(Sphere_2_Initial, T.Vector_Multiply(t_collision,Vq));
				
				//Find L = P(t_collision) - Q(collision) / |P(t_collision) - Q(collision)|
				//
				// Ensure Vector L is pointing towards Q from P as C is calculated using P(t_collision) [this is just choice]
				//
				double[] L = T.Unit_Vector(T.Vector_Subtract(Q,P));
								
				//Initialise Collision_Information. 
				//
				// Only 1 point of collision can occur for sphere to sphere collisions
				Collision_Information C = new Collision_Information(1);
				
				// Calculate Collision Point
				//
				// C = P(t_collision) + R1*L
				C.Point[0] = T.Vector_Add(P,T.Vector_Multiply(R1, L));
				C.time = t_collision;
				
				T.Print_Vector(("t_collision: " + t_collision + " ,  Collision Point"),C.Point[0]);
				
				return C;
			}
				
			//Calculates Collision point between point (Sphere r = 0) and 1D Plane - Moving or Stationary
							
			static Collision_Information Point_to_Point(double[] Sphere_1_Initial, double[] Sphere_1_Velocity, double[] Sphere_2_Initial, double[] Sphere_2_Velocity)
			{
				return Sphere_to_Sphere(Sphere_1_Initial, Sphere_1_Velocity, Sphere_2_Initial, Sphere_2_Velocity, 0, 0);
			}
			
			//Calculates Collision point between point (Sphere r = 0) and 1D Plane - Moving or Stationary
			
			static Collision_Information Plane_1D_to_1D_Plane(double[] Sphere_1_Initial, double[] Sphere_1_Velocity, double[] Sphere_2_Initial, double[] Sphere_2_Velocity)
			{
				return Sphere_to_Sphere(Sphere_1_Initial, Sphere_1_Velocity, Sphere_2_Initial, Sphere_2_Velocity, 0, 0);
			}
			
			//Calculates Collision point between Sphere and 1D Plane - Moving or Stationary
							
			static Collision_Information Sphere_to_Plane_1D(double[] Sphere_1_Initial, double[] Sphere_1_Velocity, double[] Sphere_2_Initial, double[] Sphere_2_Velocity, double R1)
			{
				return Sphere_to_Sphere(Sphere_1_Initial, Sphere_1_Velocity, Sphere_2_Initial, Sphere_2_Velocity, R1, 0);
			}
			
			//Calculates Collision point between Sphere and Point - Moving or Stationary
					
			static Collision_Information Sphere_to_Point(double[] Sphere_1_Initial, double[] Sphere_1_Velocity, double[] Sphere_2_Initial, double[] Sphere_2_Velocity, double R1)
			{
				return Sphere_to_Sphere(Sphere_1_Initial, Sphere_1_Velocity, Sphere_2_Initial, Sphere_2_Velocity, R1, 0);
			}
			
}
