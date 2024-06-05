/*SAS CODE*/
%macro LR (data=, /* input data set */
	var1=, /* continuous covariates in the logistic regression */
	var2=, /* categorical covariates in the logistic regression */
	p1=, /* number of continuous covariates in the logistic regression */
	p2=, /* number of categorical covariates in the logistic regression */
	resp=, /* binary response variable in the logistic regression */
	ntrt=); /* position of the  treatment variable in the categorical covariates */
	
	
	
	/* delete the observations with missing values in either numeric or categorical covariates */

	data newdata;
		set &data;

		if cmiss(of &var1 &var2) > 0 then delete;
	run;

	/* extract the design matrix and variance-covariance matrix of the regression and store them in data1 and parms respectively*/
	/* cases when there is no continuous covariates */

	%if &p1 = 0 %then
		%do;

			proc transreg data = newdata dummy noprint;
				model identity(&resp) = class(&var2/zero = first);
				output out = data1(drop= _: &var2 &resp) replace;
			run;

			proc logistic data = newdata desc covout outest=parms(drop=_:) noprint;
				class &var2 / param=ref ref=first;
				model &resp = &var2;
			run;

		%end;

	/* cases when both categorical and continuous covariates exist */
	%else
		%do;

			proc transreg data = newdata dummy noprint;
				model identity(&resp) = class(&var2/zero=first) identity(&var1);
				output out = data1(drop=_: &var2 &resp) replace;
			run;

			%put &var2 &var1;
			
			proc logistic data = newdata desc covout outest = parms(drop=_:) noprint;
				class &var2 / param=ref ref=first;
				model &resp  = &var2 &var1;
			run;

		%end;

	/* estimate the proportion difference and standard error to construct the 95% confidence interval*/
	proc iml;
		use parms;
		  read all into params;
		close parms;
		
		use data1;
		  read all into xmat;
		close data1;
		
		
		p = nrow(params) - 2;
		cov = params[2:(p+2), ];
		psize = nrow(xmat);
		ntrt = &ntrt + 1;
		
		xt = xmat;
		xt[,ntrt] = J(psize, 1, 1);
		axT = xt*params[1,]`;

        xc = xmat;
	    xc[,ntrt] = J(psize,1,0); 
	    axC= xc*params[1,]`;
	    
	    
		pderT = J(psize, 1, 0);
		pderC = J(psize, 1, 0);
		pderivT = J(psize, (p + 1), 0);
		pderivC = J(psize, (p + 1), 0);
		
		gt = J(1, (p + 1), 0);
		gc = J(1, (p + 1), 0);

		do k = 1 to psize;
			pderT[k] = exp(axT[k])/(1 + exp(axT[k]));
			pderC[k] = exp(axC[k])/(1 + exp(axC[k]));
			
			pderivT[k,] = pderT[k] * (1 - pderT[k]) * xt[k,];
			pderivC[k,] = pderC[k] * (1 - pderC[k]) * xc[k,];
		end;
		
		difb = sum(pderT - pderC) / psize;
		pT = sum(pderT) / psize;
		pC = sum(pderC) / psize;

		do k = 1 to (p + 1);
			gt[k] = sum(pderivT[,k])/psize;
			gc[k] = sum(pderivC[,k])/psize;
		end;
		
		se = sqrt(gt*cov*gt` + gc*cov*gc` -2*gt*cov*gc`);
		upperb = difb + 1.96*se;
		lowerb = difb - 1.96*se;
		
		store difb se pT pC lowerb upperb;
		
		xys =    difb || se || pT || pC || lowerb || upperb ;
        cname = {"diff" "se" "pt" "pC" "lower" "upper" };
		
		
		/* MB: addition to the macro to store the results as data */
		create geout from xys [ colname=cname ];
        append from xys;
		
		quit;
		
	%mend;

	/* print out the LR estimation results */
	%macro printLR (data=, var1=, var2=);
		title 'Study:' "&data";

	proc iml;
		load difb se pT pC lowerb upperb;
		
		print,  'LR estimation of study' "&data", ;
		print,  '*******************************************************', ;
		print,  'Covariates: ' "&var1" "&var2" , ;
		print , 'PropDiff = ' difb , ;
		print,  'PropT = ' pT , ;
		print , 'PropC = ' pC , ;
		print,  'SE = ' se , ;
		print,  '95% Confidence Interval= [' lowerb ',' upperb ']', ;
		print, '*********************************************************', ;
		quit;
	%mend;