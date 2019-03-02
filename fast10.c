/* fast sampling process for calculating the inclusion frequencies f_i */
/* for simple estimators for network sampling */
/* 
   Author:  Steven K. Thompson
   Copyright (C) 2019, all rights reserved
*/
	      /* idealized design 10 is within the finding sample */
	      /*  otherwise it's like idealized 1 */
	      /*** bernoulli tracing with adaptive acquisition **** 
	       *** and attrition processes ************************
	       ****************************************************/
	      if (idealized == 10)
		{
		  /* printf("t %d 0 naad %d\n", t, naad); */
		  /* save old values add up outdegree from sample */
		  aadoutdegree = 0;
		  aadtotaldegree = 0;
		  maxselects = 0.0;
		  count = 0;
		  tmp = first;
		  while(tmp != NULL)
		    {
		      /* save old value of sample membership */
		      /* here, aad[0] should always be the current
			 value, whle aad[1], aad[2], ... aad[aadlags]
			 are that many lags in the past */
		      /* even if naadlags is 0, i still need to save
			 yesterday's value as aad[1] to prevent order-dependent
			 explosions within a time step */
		      if (naadlags == 0)
			{
			  tmp->aad[1] = tmp->aad[0];
			  edgetmp = tmp->adjfirst;
			  while (edgetmp != NULL)
			    {
			      edgetmp->aadimpute[1] = edgetmp->aadimpute[0]; 
			      edgetmp = edgetmp->adj;
			    }
			}
		      else
			{
			  for (i = naadlags; i > 0;i--)
			    {
			      tmp->aad[i] = tmp->aad[i-1];
			    }
			}
		      /* add up active sample out degree if needed */
		      if (1)
			{
			  if (tmp->aad[0] == 1)
			    {
			      edgetmp = tmp->adjfirst;
			      while ( edgetmp != NULL )
				{
				  aadtotaldegree = aadtotaldegree + 1;
				  if (edgetmp->idaddress->found[0] > 0 && edgetmp->idaddress->aad[0] == 0)
				    {
				      aadoutdegree = aadoutdegree + 1;
				    }
				  if (edgetmp->idaddress->found[0] > 0 && aadreplace == 1 && edgetmp->idaddress->aad[0] == 2)
				    {
				      aadoutdegree = aadoutdegree + 1;
				    }
				  edgetmp = edgetmp->adj;
				}
			    }
			}
		      tmp = tmp->next;
		    }
		  /* printf("conductance %f out deg %d total deg %d\n", (float) aadoutdegree/aadtotaldegree, aadoutdegree, aadtotaldegree); */

		  /* temporarily sample with replacement if there are no other links to follow */
		  tempreplace = aadreplace;
		  if (aadreplace == 0 && aadoutdegree == 0 && aadfront == 1)
		    {
		      tempreplace = 1;
		    }
		  /* printf("  outdegree %d\n", aadoutdegree); */

		  /* link tracing */
		  /* printf("     1 naad %d\n", naad); */
		  /* only trace if constant rate or sample size < target */
		  if (aadfront == 0 || (aadfront == 1 && naadtarget - naad > 0))
		    {
		      /* set tracing rate */
		      if (aadfront == 0)
			/* constant tracing rate */ 
			{
			  u = paadtrace;
			}
		      else /* if aadfront == 1 */
			/* control rate to give expected target sample size */
			{
			  if (naadtarget - naad > 0)
			    {
			      if (tempreplace != aadreplace)
				{
				  u = 1.0/aadoutdegree;
				}
			      else
				{
				  u = ((float) naadtarget - naad)/aadoutdegree;
				  if (u > 1.0) u = 1.0;
				}
			    }
			  else u = 0.0;
			}

		      /* printf("IN FAST 10 aadfront %d u %f\n", aadfront, u); */

		      /* now trace using rate u */ 
		      /* printf("ntarget - n %d, outdeg %d, u %f\n", (naadtarget-naad), aadoutdegree, u); */
		      tmp = first;
		      while (tmp != NULL)
			{
			  /* use the past value aad[1] so it's not 
			     recursively traced if it turns 1 at this time */
			  if (tmp->aad[1] == 1)
			    {
			      edgetmp = tmp->adjfirst;
			      while ( edgetmp != NULL )
				{
				  /* make option for findtraced either
				     here or in finding design */
				  if (edgetmp->idaddress->found[0] > 0 && edgetmp->findtraced == 1)
				  /* if (edgetmp->idaddress->found[0] > 0) */
				    {
				      if (edgetmp->idaddress->aad[0] == 0 || (edgetmp->idaddress->aad[0] == 2 && tempreplace == 1))
					{
					  if (coupons == 0)
					    {
					      /* if ((edgetmp->safersex == 0 && runif(mts3) < u) || (edgetmp->safersex > 0 && runif(mts3) < u * safersexeffect)) */
					      if (runif(mts3) < u)
						{
						  edgetmp->idaddress->aad[0] = 1;
						  edgetmp->tracing = t;
						  naad = naad + 1;
						} 
					    }
					  if (coupons > 0)
					    {
					      if (tmp->coupon > 0 && runif(mts3) < u)
						{
						  edgetmp->idaddress->aad[0] = 1;
						  edgetmp->tracing = t;
						  edgetmp->idaddress->coupon = coupons;
						  tmp->coupon = tmp->coupon - 1;
						  naad = naad + 1;
						  /* printf("check 2  u %f, coupon %d id %d edge id %d\n",u, tmp->coupon, tmp->id, edgetmp->id); */
						}
					    }
					} /* if edgetmp idaddress aad[0] etc */
				    } /* if edgetmp->idaddress->found > 0 */


				  edgetmp = edgetmp->adj;
				}
			    }
			  tmp = tmp->next;
			}

		      /* tracing of fast design to and from imputed node */
		      if (ifindimputing == 1)
		      	{
			  tmp = first;
			  while(tmp != NULL)
			    {
			      if (tmp->found[0] == 1)
				{
				  edgetmp = tmp->adjfirst;
				  while(edgetmp != NULL)
				    {
				      if (edgetmp->findimpute > 0)
					{
					  /* trace from host to imputed */ 
					  if (tmp-> aad[1] == 1 && edgetmp->aadimpute[0] == 0)
					    {
					      if (coupons == 0)
						{
						  /* if ((edgetmp->safersex == 0 && runif(mts3) < u) || (edgetmp->safersex > 0 && runif(mts3) < u * safersexeffect)) */
						  if (runif(mts3) < u)
						    {
						      edgetmp->aadimpute[0] = 1;
						      edgetmp->tracing = t;
						      naad = naad + 1;
						    } 
						}
					      if (coupons > 0)
						{
						  if (tmp->coupon > 0 && runif(mts3) < u)
						    {
						      edgetmp->idaddress->aad[0] = 1;
						      edgetmp->tracing = t;
						      edgetmp->idaddress->coupon = coupons;
						      tmp->coupon = tmp->coupon - 1;
						      naad = naad + 1;
						      /* printf("check 2  u %f, coupon %d id %d edge id %d\n",u, tmp->coupon, tmp->id, edgetmp->id); */
						    }
						}
					    } /* if edgetmp idaddress aad[0] etc */
					  /* trace from imputed to host */ 
					  else if (tmp-> aad[1] == 0 && edgetmp->aadimpute[0] > 0)
					    {
					      if (runif(mts3) < u)
						{
						  tmp->aad[0] = 1;
						  tmp->coupon = coupons;
						  edgetmp->tracing = t;
						  naad = naad + 1;
						} 
					    }
					} /* if edgetmp->findimpute > 0 */
				      edgetmp = edgetmp->adj;		    
				    }
				} /* if tmp-> found == 1 */
			      tmp = tmp->next;
			    }
			} /* if ifindimputing */

		      /* end tracing of fast design to and from imputed node */
		    } /* if aadfront = 0 or (1 and less than target) */
		  /* attrition from aad sample */
		  /* printf("     2 naad %d\n", naad); */
		  /* only do if constant attrition rate or sample size > target */
		  if (aadback == 0 || (aadback == 1 && naad - naadtarget > 0))
		    {
		      /* set attrition rate */
		      if (aadback == 0)
			/* constant attrition rate */
			u = paadremove;
		      else /* if aadback == 1 */ 
			/* attrition rate to give expected target sample size */
			{
			  if (naad - naadtarget > 0)
			    {
			      u = ((float) naad - naadtarget)/naad;
			    }
			  else u = 0.0; 
			}
		      tmp = first;
		      while(tmp != NULL)
			{
			  if (coupons == 0) /* no coupon program */
			    {
			      if (tmp->aad[0] == 1 && runif(mts3) < u)
				{
				  tmp->aad[0] = 2;
				  /* i may need to make transmission -30 here */
				  /* or use aad[0] to not delete on find step */
				  naad = naad - 1;
				  naad2 = naad2 + 1;
				}
			    }
			  else /* coupon program */
			    {
			      /* remove from active sample if out of coupons */
			      /* or by attrition rate */ 
			      if (tmp->aad[0] == 1 && (tmp->coupon == 0 || runif(mts3) < u))
				{
				  tmp->aad[0] = 2;
				  naad = naad - 1;
				  naad2 = naad2 + 1;
				  if (tmp->coupon < 0)
				    printf("ERROR fewer than 0 coupons\n");
      				}
			    }
			  /* attrition of imputed node */
			  edgetmp = tmp->adjfirst;
			  while(edgetmp != NULL)
			    {
			      if (edgetmp->aadimpute[0] == 1 && runif(mts3) < u)
				{
				  edgetmp->aadimpute[0] = 0;
				}
			      edgetmp = edgetmp->adj;
			    }
			  tmp = tmp->next;
			}
		    } /* if aadback == 0 etc */
		  /* printf("     3 naad %d\n", naad); */


		  /* random selections, reseeds */
		  tmp = first;
		  while(tmp != NULL)
		    {
		      /*  random reseeds */
		      if (aadreseeddesign == 1)
			{
			  aadreseedprob = paadrandom;
			} /* random reseeds */
		      /* reseeds with probability proportional to degree */
		      else if (aadreseeddesign == 2)
			{
			  /* want average of pp degree to be pfindingreseed */
			  /* c * sum_nodes tmp->degree/N_t = pfindreseed */
			  /* c * 2.0 * nedges/Nt = pfindreseed */
			  /* c = 0.5 * pfindreseed * Nt/nedges */
			  aadreseedprob = 0.5 * paadrandom * Nt * tmp->degree/nedge;
			}
		      /* select reseeds from nodes not in current sample */
		      if (tmp->aad[0] != 1 && tmp->found[0] > 0)
			{
			  /* the factor 10.0 because this design needs higher
			     reseed prob than a without replacement one */
			  if (runif(mts3) < 1.0 * aadreseedprob)
			    {
			      tmp->aad[0] = 1;
			      /* tmp->aadselects[0] = tmp->aadselects[0] + 1; */
			      naad = naad + 1;
			      if (coupons > 0)
				{
				  tmp->coupon = coupons;
				}
			      /* if (coupons > tmp->degree) */
			      /* 	{ */
			      /* 	  tmp->coupon = tmp->degree; */
			      /* 	} */
			      /* else */
			      /* 	{ */
			      /* 	  tmp->coupon = coupons; */
			      /* 	} */
			    }
			}
		      tmp = tmp->next;
		    }

		  /* averaging aad values */
		  tmp = first;
		  while(tmp != NULL)
		    {
		      /* moving average over nlags + 1 steps */


		      if (aadar == 0)
			{
			  /* average aad over naadlags + 1 steps */
			  tmp -> aadlagav = 0.0;
			  for (i=0; i <= naadlags; i++)
			    {
			      /* printf("%d ", tmp->aad[i]); */
			      if (tmp->aad[i] == 1)
				{
				  tmp->aadlagav = tmp->aadlagav + 1.0;
				}
			    }
			  tmp->aadlagav = tmp->aadlagav/(1.0 + naadlags);
			}  /* if aadar == 0 */
		      else if (aadar == 1)
			/* autoregress */ 
			{
			  /* tmp->aadlagav = tmp->aadlagav/(1.0 - aadarphi); */
			  if (tmp->aad[0] == 1)
			    {
			      tmp->aadlagav = tmp->aadlagav * aadarphi + (1.0 - aadarphi) * tmp->aad[0];
			    }
			  else /* aad[0] = 0 or 2 both treated as 0 */
			    {
			      tmp->aadlagav = tmp->aadlagav * aadarphi;
			    }
			  /* tmp->aadlagav = tmp->aadlagav * (1.0 - aadarphi); */
			}
		      /* cumulative mean of aad */
		      else if (aadar == 2)
			/* autoregress */ 
			{
			  /* tmp->aadlagav = tmp->aadlagav/(1.0 - aadarphi); */
			  if (tmp->aad[0] == 1)
			    {
			      tmp->aadlagav = (-1.0 + t)/t * tmp->aadlagav + 1.0/t * tmp->aad[0];
			      /* make aadselects[0] 1 if aad[0] is 1 */
			      /* for use elsewhere in simple estimation */
			      /* tmp->aadselects[0] = 1; */
			    }
			  else /* aad[0] = 0 or 2 both treated as 0 */
			    {
			      tmp->aadlagav = tmp->aadlagav * (-1.0 + t)/t;
			      /* tmp->aadselects[0] = 0; */
			    }
			  /* tmp->aadlagav = tmp->aadlagav * (1.0 - aadarphi); */
			}
		      if (tmp->aadlagav > maxselects)
			{
			  /* printf("aadlagav %f maxselects %f id %d born %d\n", tmp->aadlagav, maxselects, tmp->id, tmp->eventtime[0]); */
			  maxselects = tmp->aadlagav;
			}
		      /* if (tmp->found > 0) */
		      /* 	{ */
		      /* 	  printf("aadlagav in design 10 %f found %d\n", tmp->aadlagav/maxselects, tmp->found); */
		      /* 	} */

		      /* put aadlagav in a vector for sorting */
		      if (iflamesort == 1)
			{
			  flamesort[count] = tmp->aadlagav;
			  count = count + 1;
			}
		      if (iflamesort == 2 && tmp->bug == 0)
			{
			  flamesort[count] = tmp->aadlagav;
			  count = count + 1;
			}
		      if (iflamesort == 3 && tmp->bug == !0)
			{
			  flamesort[count] = tmp->aadlagav;
			  count = count + 1;
			}
		      /* rank (descending) below which we treat */
		      flamemintreat = (int) fmin(count, 100);
		      /* if (t>0) */
		      /* printf("naadlagav %f %f\n", tmp->aadlagav, tmp->aadlagav * (1.0 - aadarphi));  */
		      tmp = tmp->next;
		    }
		  if (iflamesort)
		    {
		      /* printf("count %d Nt %d nobug %d nbug %d\n", count, Nt, Nt-nbug, nbug); */
		      qsort(flamesort, count, sizeof(double), compare);
		    }

		  /* printf("maxselects %f\n", maxselects); */
		  /* printf("FLAMESORT\n"); */
		  /* for (j=0; j < 20; j++) */
		  /*   { */
		  /*     printf("%f ", flamesort[j]); */
		  /*   } */
		  /* printf("\n"); */
		  /* maxselects = 1.0; */
		} /* if idealized == 10 */
	      /************************************/
