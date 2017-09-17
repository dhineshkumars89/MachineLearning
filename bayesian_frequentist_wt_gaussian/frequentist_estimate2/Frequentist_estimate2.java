import java.util.Random;
import java.util.concurrent.TimeUnit;

public class Frequentist_estimate2 {

	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		
		 int a = 0,b = 0,c = 0,d = 0,e = 0;
		 for(int j=0;j < 10000; j++)
		 {
			 Random rnd = new Random();
			 //rnd.setSeed(2000);
			
			 float val;
			 //String str = "";
			 StringBuilder str = new StringBuilder(3100);
			 float numa = 0;
			 for(int i=0;i < 3100; i++)
			 {
				 val = rnd.nextFloat();
				 if(val <= 0.1)
				 {
					 str.append('a');
					 numa++;
				 }
				 else
				 {
					 str.append('b');
				 }
			 }
			 
			 float prob = numa/3100;
			 if( prob < 0.08)
				 a++;
			 else if(prob < 0.09)
				 b++;
			 else if((prob >= 0.09) && (prob <= 0.11))
				 c++;
			 else if(prob > 0.12)
				 e++;
			 else if(prob > 0.11)
				d++;
			 else
			 {
				 //do nothing
			 }
		 }
		 //System.out.println(str);
		 System.out.println("In "+a+" of the simulations p(c = 'a') < 0.08.");
		 System.out.println("In "+b+" of the simulations p(c = 'a') < 0.09.");
		 System.out.println("In "+c+" of the simulations p(c = 'a') is in the interval [0.09, 0.11].");
		 System.out.println("In "+d+" of the simulations p(c = 'a') > 0.11.");
		 System.out.println("In "+e+" of the simulations p(c = 'a') > 0.12.");
		 
		 		 
		 long endTime   = System.currentTimeMillis();
		 long totalTime = endTime - startTime;
		 //System.out.println(totalTime);
		 //System.out.println(TimeUnit.MILLISECONDS.toSeconds(totalTime));
	}

}
