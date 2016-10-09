
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;

//Rotate Body
//
// 5 - Change Rotation Matrix to General Transform Matrix
//

public class Main {

	public static void main(String[] args)
	{	
		
		double[] T = {1,1,1};
		
		double[][] Box = { {0,0,4},
				{1,0,4},
				{0,1,4},
				{0.5,2,4},
				{1,1,4},
		};
		
		Body object = new Body(Box);
		
		Body.CalculateCoM();
		Body.Mass = 5;
		Body.CalculateCoM();
		Maths.Print_Vector("CoM",object.CoM);
		
		
	/*	double[] a = {1,0,0,0};
		double[] b = new double[4];
		double[] euler = {T.Degrees_to_Radians(90),0,0};
		double[] p = {0,0,1};

		double[] M = T.Euler_to_Quaternion(euler);
		T.Print_Vector("Euler", euler);
		T.Print_4D_Vector("Euler.Quaternion", M);
		double[] AxisAngle = T.Quaternion_to_AxisAngle(M);
		T.Print_4D_Vector("Euler.AxisAngle", AxisAngle);
		M = T.Quaternion_to_Euler(M);
		T.Print_Vector("Euler", M);*/

	}
	
}
