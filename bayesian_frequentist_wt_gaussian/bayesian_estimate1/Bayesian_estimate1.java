public class Bayesian_estimate1 {

	public static void main(String[] args) {
		
		String data;
		if(args.length == 0)
			data = "";
		else
			data = args[0];

		double[] theta = {0.1, 0.3, 0.5, 0.7, 0.9};
		double[] probtheta = {0.9, 0.04, 0.03, 0.02, 0.01};
		float numa = 0;
		float numb = 0;
		float prob_a_data = 0;
		//process data
		for(int i=0;i<data.length();i++)
		{
			if(data.charAt(i) == 'a')
				numa++;
			//if(data.charAt(i) == 'b')
			else
				numb++;
		}
		double[] p_data_gvn_theta = new double[theta.length];
		double posterior_prob,totalval=0;
		double total=0;
		for(int i=0;i<theta.length; i++)
		{
			 double val = ((Math.pow(theta[i],numa)) * (Math.pow((1-theta[i]),numb)));
			 p_data_gvn_theta[i] = val * probtheta[i];
			 total = total + p_data_gvn_theta[i];
		}
		for(int j=0;j<theta.length;j++)
		{
			posterior_prob = (p_data_gvn_theta[j]/total);
			System.out.println(String.format("p(m = %.1f | data) = %.4f",theta[j],posterior_prob));
			totalval = totalval + (posterior_prob * theta[j]);
			if(j== theta.length-1)
			{
				System.out.println(String.format("p(c = 'a' | data) = %.4f",totalval));
			}
		}
		
	}

}
