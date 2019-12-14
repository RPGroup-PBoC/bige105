%In this script, we will see how the poisson distribution changes with
%different number of trials and different probabilities of success. As a
%reminder, the Poisson distribution is defined as
% 
%     P(n,N) = (lambda^n / n!) * exp(-lambda)
% 		lambda = N * p
%
%where N is the number of trials, n is the number of successes, and p is the
%probability of a success with any given coin flip. For our exercise, we will
%vary N and p.

%Start out by defining some parameters. We will start by making a vector of N's
%for each plot. 
N = logspace(1, 10, 20); %10 values from 10 to 10^10
n = 0:1:100; %100 values from 0 to 100.
p = 1E-5; %Probability of a success. 

%Since we will change the value of lambda at each step, we will have to loop
%through each N value we have defined above. 
figure(1)
for i=1:length(N)
	%Now define lambda. 
	lam = N(i) * p;
	
	%Evaluate the Poisson. Since we are dealing with vectors, we must
	%divide with a ./ instead of simply /  
	poiss = ((lam.^n) ./ factorial(n)) * exp(-lam);

	%Now plot it and give it a label. 
	plot(n, poiss, '-', 'DisplayName', ['lambda = ' num2str(lam)])
	hold on;
end

%Now we just have to add the most important parts -- labels and a legend. 
xlabel('number of successes')
ylabel('probability')
legend('-DynamicLegend') %This will apply the labels we defined above. 
hold off;

%We can see that as lamda becomes large, the poisson becomes more reminiscient
%of a gaussian distribution. Conversely, when lambda becomes very small, the
%distribution looks much more like a binomial distribution.
