
//Contains information regarding body in world

public class Body {

	//Body Name
	static String Name = "";
	
	//Mechanical Information - Rotation is Sequenced as Quaternion
	static double[] Angular_Velocity = new double[4];
	
	//Linear Velocity for current Animation Frame
	static double[] Linear_Velocity = new double[4];
	
	//Acceleration of Body in Motion
	public double Linear_Acceleration;
	public double Angular_Acceleration;
	
	static double Spring_Constant;
	
	//Mass for Non-Articulated Bodies is assumed uniform
	static double Mass;
	static double[] CoM;
	
	//Geographical Information
	static double[][] Body_Vertices;
	static double[][] Articulated_Points;
	
	//Lighting Information
	static double[] Lighting_Properties;
	static double[] Material_Properties;
	
	Body(double[][] Body)
	{
		//Initialise quaternion rotation to [1 0 0 0]
		Angular_Velocity[0] = 1;
		Angular_Velocity[1] = Angular_Velocity[2] = Angular_Velocity[3] = 0;
		
		//Default Linear Vector of CoM to 0
		Linear_Velocity[0] = Linear_Velocity[1] = Linear_Velocity[2] = Linear_Velocity[3] = 0;
		
		//Store Body Vertice Map
		Body_Vertices = Body;
		
		CoM = null;
		Mass = 0;
		
	}
	
	//Calculate bodies centre of mass
	static void CalculateCoM()
	{
		//CoM = Summation of (vertex*mass @ vertex)/sum of all vertex masses
		
		if(Mass == 0)
		{
			System.out.println("Must know Body Mass before calculating CoM");
			return;
		}
			
		double[] numerator = {0,0,0};
		double denominator = 0;
		
		int i;
		for(i = 0; i < Body_Vertices.length; i++)
		{
			numerator = Maths.Vector_Add(numerator,Maths.Vector_Multiply(Mass,Body_Vertices[i]));
			denominator += Mass;
		}
		
		CoM = Maths.Vector_Multiply( (1/denominator), numerator);
		
	}

	//Change Name of body
	//
	static void BodyName(String s)
	{
		Name = s;
	}
	
	//
	// Translates non-articulated body by CoM Linear Vector
	// Does not use Affine Transformation Matrix
	//
	static void TranslateBody()
	{
		Maths T = new Maths();

		int i;
		for(i = 0; i < Body_Vertices.length; i++)
		{
			Body_Vertices[i] = T.Vector_Add(Body_Vertices[i], Linear_Velocity);
		}	
	}
	
	//
	// Convert Quaternion Rotation Sequence to AxisAngle.
	// Use AxisAngle & Rodriguez formulae to apply rotation
	//
	
	static void RotateBody_AxisAngle()
	{
		// Angle of AxisAngle in Radians
		//
		
		Maths T = new Maths();
		
		//Create AxisAngle from rotation sequence
		double[] AxisAngle = T.Quaternion_to_AxisAngle(Angular_Velocity);
		
		//Apply Rotation to each body vertex
		int i;
		for(i = 0; i < Body_Vertices.length; i++)
		{
			Body_Vertices[i] = T.AxisAngle_Rotation(AxisAngle, Body_Vertices[i]);
		}
		
	}
	
	//
	// Convert Quaternion Rotation Sequence to Quaternion Rotation Matrix
	// Use Quaternion Matrix to apply rotation to body
	//
	
	public void RotateBody_Quaternion()
	{
		Maths T = new Maths();
		
		//Store Body Vertice Map
		double[][] points = new double[Body_Vertices.length][4];
		points = Body_Vertices;
				
		//Set 4th dimension of Body Vertices to enable Translation in affine transform
		for(int i = 0; i < Body_Vertices.length; i++)
		{
			points[Body_Vertices.length][4] = 1;
		}
		
		//Create quaternion Rotation Matrix
		double[][] R = T.Quaternion_to_QuaternionMatrix(Angular_Velocity);
		
		//Rotate Body
		double[][] temp = new double[Body_Vertices.length][3];
		temp = T.Matrix_Multiply(R,points);
		
		Body_Vertices = temp;
	}
	
	//
	// Convert Quaternion Rotation Sequence to Euler Matrix.
	// Use Euler Matrix to apply rotation to body.
	//
	public void RotateBody_Euler()
	{
		Maths T = new Maths();
		
		//Store Body Vertice Map
		double[][] points = new double[Body_Vertices.length][4];
		points = Body_Vertices;
				
		//Set 4th dimension of Body Vertices to enable Translation in affine transform
		for(int i = 0; i < Body_Vertices.length; i++)
		{
			points[Body_Vertices.length][4] = 1;
		}
		
		//Create Euler Matrix
		double[][] R = T.Euler_to_EulerMatrix(T.Quaternion_to_Euler(Angular_Velocity));
		
		double[][] temp = new double[Body_Vertices.length][3];
		temp = T.Matrix_Multiply(R,points);
		
		Body_Vertices = temp;
		
	}
	
	// Set Linear Velocity for [X Y Z]
	public void SetLinearVelocity(double[] V)
	{
		Linear_Velocity = V;
		Linear_Velocity[3] = 0;
	}
	
	// rotation angle = Degrees
	//
	// Multiplies Current Rotation Quaternion by next Quaternion in sequence (like glRotatef)
	//
	public void AddRotationToSequence(double rot, double a, double b, double c)
	{
		Maths T = new Maths();
		double[] AxisAngle = {rot, a, b, c};
		
		rot = T.Degrees_to_Radians(rot);
		AxisAngle = T.AxisAngle_to_Quaternion(AxisAngle);
		Angular_Velocity = T.Quaternion_Multiplication(AxisAngle,Angular_Velocity);
		
	}
	
	//rotation angle = Degrees
	//
	//Resets Rotation Quaternion and Creates New Rotation Sequence	
	//
	public void SetRotation(double rot, double a, double b, double c)
	{
		Maths T = new Maths();
		double[] AxisAngle = {T.Degrees_to_Radians(rot), a, b, c};
		
		Angular_Velocity = T.AxisAngle_to_Quaternion(AxisAngle);
	}
	
}
