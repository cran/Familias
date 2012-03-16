/* $Id: alsys.C,v 1.11 2001/04/12 16:00:35 mostad Exp $ */

#include <cstring>
#include <iostream>

#include "family.h"
#include "odds.h"
#include "alsys.h" 
#include "table.h"
#include "person.h"
#include "special.h"

///////////////////////////////////////////
//Relating to dataitem:
///////////////////////////////////////////

void dataitem::remove_next(dataitem* dit) {
    if (!next) return; 
    if (dit==next) next = next->next;
    else next->remove_next(dit);
}

///////////////////////////////////////////
//Relating to allelesystem:
///////////////////////////////////////////


allelesystem::allelesystem(char* sname, 
			   double rateFemale, 
			   double rateMale, 
			   int mutationModFemale, 
			   int mutationModMale, 
			   int n_poss, 
			   double mutRangeFemale, 
			   double mutRangeMale) :
   Systemname(sname), 
   mutationrateFemale(rateFemale), 
   mutationrateMale(rateMale), 
   mutationModelFemale(mutationModFemale), 
   mutationModelMale(mutationModMale), 
   n_possibilities(n_poss), 
   mutationRangeFemale(mutRangeFemale), 
   mutationRangeMale(mutRangeMale), 
   kinship(0), 
   next(0),
   n_alleles(0), 
   name(0), 
   probability(0), 
   hasSilentAllele(0),
   silentAllele(0),
   data(0), 
   result(1),
   n_dataalleles(0), 
   index(0), 
   dataprobability(0), 
   dataprobmatrixFemale(0),
   dataprobmatrixMale(0),
   recalc_data(1) {}

allelesystem::~allelesystem() {
    int i;
    delete[] Systemname;
    for (i=0; i<n_alleles; i++) 
	delete[] name[i];
    delete[] name;
    delete[] probability;
    delete[] index;
    if (data) delete data;

    if (n_dataalleles) {
	delete[] dataprobability;
	for (i=0; i<n_dataalleles; i++) 
	{
	    delete[] dataprobmatrixFemale[i];
	    delete[] dataprobmatrixMale[i];
	}
	delete[] dataprobmatrixFemale;
	delete[] dataprobmatrixMale;
    }
}

void allelesystem::remove_next(allelesystem* s) {
    if (!next) return; 
    if (s==next) next = next->next;
    else next->remove_next(s);
}

int allelesystem::new_mutationrate(double rateFemale,
				   double rateMale, 
				   int mutationModFemale, 
				   int mutationModMale, 
				   int n_poss, 
				   double mutRangeFemale, 
				   double mutRangeMale, 
				   int info, 
				   int& error) 
{
   if (rateFemale == mutationrateFemale && 
       rateMale   == mutationrateMale   &&
       mutationModelFemale == mutationModFemale &&
       mutationModelMale   == mutationModMale   &&
       //(mutationModelFemale > 1 || n_poss == n_possibilities) &&
       //(mutationModelMale   > 1 || n_poss == n_possibilities) &&
       (mutationModelFemale < 2 || mutRangeFemale == mutationRangeFemale) &&
       (mutationModelMale   < 2 || mutRangeMale   == mutationRangeMale))
   {
      return 0;
   }

   
   
// THESE TESTS SHOULD NOT ME NECESSARY NOW:    
//   if ((mutationModelFemale != 2 || mutationModelMale != 2) 
//       && n_poss < n_alleles) 
//   {
//      if (info>0)
//	 cout<<"ERROR: There are "<<n_alleles
//	     <<" registered alleles in the system \""
//	     <<Systemname
//	     <<"\".\nThe number of possible alleles must be "
//	     <<"at least that number.\n";
//      error = 1; 
//      return 0;
//   }
//   if ((mutationModelFemale != 2 || mutationModelMale != 2)
//       && n_poss == n_alleles) 
//   { //The sum of the frequencies in the system must be 1.
//      double sum = 0;
//      for (int i = 0; i<n_alleles; i++) sum += probability[i];
//      if (sum <= 0.999 || sum >= 1.001) 
//      {
//	 if (info>0)
//	    cout<<"ERROR: The number of possible alleles cannot "
//		<<"be set equal to the\n"
//		<<"actual number of alleles unless "
//		<<"their frequencies sum to 1.\n";
//	 error = 1; 
//	 return 0;
//      }
//   }

   mutationrateFemale    = rateFemale;
   mutationrateMale      = rateMale;
   mutationModelFemale   = mutationModFemale; 
   mutationModelMale     = mutationModMale; 
   n_possibilities = n_poss;
   mutationRangeFemale   = mutRangeFemale; 
   mutationRangeMale     = mutRangeMale; 
   recalc_data = 1;
// Removed 2012-03-08
//   if (info>1)
//   {
//      cout<<"The system "<<Systemname<<" now has female mutation rate "
//	  <<rateFemale<<",\nmale mutation rate"<<rateMale
//	  <<",\nand "<<n_poss<<"  possible alleles.\n";
//   }
   return 1;
}

void allelesystem::setMutation(int mutMod, 
			       double mutRange)
{
   mutationModelFemale = mutMod; 
   mutationModelMale = mutMod; 
   mutationRangeFemale = mutRange; 
   mutationRangeMale = mutRange; 
   recalc_data = 1; 
}

void allelesystem::setKinship(double kship)
{
   kinship = kship; 
   recalc_data = 1; 
}

int allelesystem::get_number_of_systems()
{
   if (!next) return 1; 
   return next->get_number_of_systems() + 1; 
}

int allelesystem::add_allele(char* allelename, 
			     double prob, 
			     int info, 
			     int& error) 
{
   int i;
   if (prob<=0) {
// Removed 2012-03-08
//      if (info>0)
//	 cout<<"ERROR: The given probability "<<prob
//	     <<" must be greater than 0!\n";
      error = 1; 
      delete[] allelename;
	  return 0;
   }
   double sum = 0;
   for (i=0; i<n_alleles; i++) sum += probability[i];
   for (i=0; i<n_alleles; i++)
      if (strcmp(name[i], allelename)==0) {
	 if (prob==probability[i]) 
	 {
	    delete[] allelename;
		return 0;
	 }

// THIS TEST CREATES PROBLEMS WHEN USED TOGETHER WITH THE 
// EditAlleleSystem command of FamInterface, SO IT IS REMOVED: 
//
//	 if ((mutationModelFemale != 2 || mutationModelMale != 2) 
//	     && n_alleles==n_possibilities) 
//	 {
//	    if (info>0)
//	       cout<<"ERROR: Resetting the probability of the allele \""
//		   <<allelename<<"\" in the system \""<<Systemname<<"\"\n"
//		   <<"is impossible as long as the number of alleles in the "
//		   <<"system\nis equal to the possible number of alleles.\n";
//	    error = 1; 
//	    delete[] allelename;
//		return 0;
//	 }


	 if (sum + prob - probability[i] >= 1.001) 
	 {
// Removed 2012-03-08
//	    if (info>0)
//	       cout<<"ERROR: Setting the probability "<<prob
//		   <<" to the allele "<<allelename
//		   <<" gives a\nprobability sum greater than 1 in "
//		   <<"the system \""<<Systemname<<"\".\n";
	    error = 1; 
	    delete[] allelename;
		return 0;
	 }
	 probability[i] = prob;
	 recalc_data = 1;
// Removed 2012-03-08
//	 if (info>1)
//	    cout<<"The allele \""<<allelename<<"\" in the system \""
//		<<Systemname<<"\" gets the probability "<<prob<<".\n";
	 delete[] allelename;
	 return 1;
      }
//   if ((mutationModelFemale < 2 || mutationModelMale < 2)
//       && n_alleles==n_possibilities) 
//   {
//      if (info>0)
//	 cout<<"ERROR: Before another allele can be added "
//	     <<"to the system "<<Systemname<<",\n"
//	     <<"the total number of possible alleles in "
//	     <<"the system must be increased from "
//	     <<n_possibilities<<".\n";
//	  error = 1; 
//      delete[] allelename;
//      return 0;
//   }
   if (sum + prob > 1.001) 
   {
// Removed 2012-03-08
//      if (info>0)
//	 cout<<"ERROR: Setting the probability "<<prob<<" to the allele "
//	     <<allelename<<" gives a\nprobability sum greater than 1 in "
//	     <<"the system \""<<Systemname<<"\".\n";
      error = 1; 
      delete[] allelename;
	  return 0;
   }

//When the program is run from the VB coat, this is unneccessary: 
//   if ((mutationModelFemale != 2 || mutationModelMale != 2) && 
//       n_alleles==n_possibilities-1 && 
//       sum + prob < 0.999) 
//   {
//      if (info>0)
//	 cout<<"ERROR: The allele \""<<allelename
//	     <<"\" is the last possible allele in\n"
//	     <<"the system \""<<Systemname
//	     <<"\", so its frequency must be "<<1-sum<<",\n"
//	     <<"and not "<<prob<<", to make the sum of the frequencies 1.\n";
//      error = 1; 
//      delete[] allelename;
//	  return 0;
//   }

   n_alleles++;
   double* newprobability = new double[n_alleles];
   char** newname = new char*[n_alleles];
   for (i=0; i<n_alleles-1; i++) 
   {
      newprobability[i] = probability[i];
      newname[i] = name[i];
   }

   delete[] probability;
   probability = newprobability;
   delete[] name;
   name = newname;
   probability[n_alleles-1] = prob;
   name[n_alleles-1] = allelename;
// Removed 2012-03-08
//   if (info>1)
//      cout<<"The allele \""<<allelename<<"\" with probability "<<prob
//	  <<" added to the system \""
//	  <<Systemname<<"\".\n";
   recalc_data = 1;
   return 1;
}

int allelesystem::add_data(person& p, char* allele1, char* allele2, 
			   int info, int& error) 
{
   int i = 0;
   for (;;) {
      if (i==n_alleles) 
      {
// Removed 2012-03-08
//	 if (info>0)
//	    cout<<"ERROR: The allele \""<<allele1<<"\" could not be found "
//	       "in the system \""<<Systemname<<"\".\n";
	 error = 1; 
	 delete[] allele1; 
	 delete[] allele2;
	 return 0;
      }
      if (strcmp(name[i], allele1)==0) break;
      i++;
   }
   int j = 0;
   for (;;) {
      if (j==n_alleles) 
      {
// Removed 2012-03-08
//	 if (info>0)
//	    cout<<"ERROR: The allele \""<<allele2<<"\" could not be found "
//	       "in the system \""<<Systemname<<"\".\n";
	 error = 1; 
	 delete[] allele1; 
	 delete[] allele2;
	 return 0;
      }
      if (strcmp(name[j],allele2)==0) break;
      j++;
   }
   dataitem* dit = data;
   while (dit) 
   {
      if (dit->p==&p) 
      {
	 if ((dit->allele1==i && dit->allele2==j) ||
		(dit->allele1==j && dit->allele2==i)) 
	 {
	    delete[] allele1; 
	    delete[] allele2;
	    return 0;
	 }
	 dit->allele1 = i;
	 dit->allele2 = j;
// Removed 2012-03-08
//	 if (info>1) 
//	    cout<<"The person \""<<p.name()<<"\" now has allele data \""
//		<<allele1<<"\" and \""<<allele2<<"\".\n";
	 recalc_data = 1;
	 delete[] allele1; 
	 delete[] allele2;
	 return 1;
      }
      dit = dit->next;
   }
   dit = new dataitem(p,i,j);
   if (data) data->append(dit);
   else data = dit;
// Removed 2012-03-08
//   if (info>1)
//      cout<<"Alleles \""<<allele1<<"\" and \""<<allele2<<"\" for \""
//	  <<p.name()<<"\" have been added to the data.\n";
   recalc_data = 1;
   delete[] allele1; 
   delete[] allele2;
   return 1;
}

int allelesystem::set_allele_as_silent(char* allelename, int& error) 
{
   int i;
   for (i=0; i<n_alleles; i++)
      if (strcmp(name[i], allelename)==0) {
		  hasSilentAllele = 1; 
		  silentAllele = i; 
		  recalc_data = 1;
		  error = 0; 
		  return 1; 
	  }
   error = 1; 
   return 0; 
}

int allelesystem::remove_allele(char* allelename, int info, int& error) 
{
   int i,j;
   for (i=0; i<n_alleles; i++)
      if (strcmp(name[i], allelename)==0) {
	 //Removing is illegal if the allele still appears in the data:
	 dataitem* dit = data;
	 while (dit) {
	    if (dit->allele1==i || dit->allele2==i) {
// Removed 2012-03-08
//	       if (info>0)
//		  cout<<"ERROR: The allele \""<<allelename<<"\" may not be "
//		      <<"removed from system \""<<Systemname<<"\".\n\""
//		      <<dit->p->name()
//		      <<"\" still has this allele as data.\n";
	       error = 1; 
	       delete[] allelename;
	       return 0;
	    }
	    dit = dit->next;
	 }
	 
	 //THIS PART WAS ADDED 21 SEPTEMBER 2009, correcting a longtime bug!!
	 dit = data; 
     while (dit) {
		 if (dit->allele1>i) dit->allele1--; 
		 if (dit->allele2>i) dit->allele2--; 
         dit = dit->next;
     }
	 //END OF ADDITION

	 delete[] name[i];
	 for (j=i+1; j<n_alleles; j++) {
	    name[j-1] = name[j];
	    probability[j-1] = probability[j];
	 }
	 if (hasSilentAllele) {
		 if (silentAllele == i) 
			 hasSilentAllele = 0;
		 else if (silentAllele > i)
			 silentAllele--; 
	 }
	 n_alleles--;
// Removed 2012-03-08
//	 if (info>1)
//	    cout<<"The allele \""<<allelename<<"\" is removed from the"
//		<<"system \""<<Systemname<<"\".\n";
	 recalc_data = 1;
	 delete[] allelename;
	 return 1;
      }
// Removed 2012-03-08
//   if (info>0)
//      cout<<"ERROR: No allele \""<<allelename<<"\" found in system \""
//	  <<Systemname<<".\n";
   error = 1; 
   delete[] allelename;
   return 0;
}

int allelesystem::remove_data(person& p, int info, int&) 
{
   dataitem* dit = data;
   while (dit) {
      if (dit->p==&p) {
	 if (dit==data) data = data->next;
	 else data->remove_next(dit);
	 dit->next = 0;
	 delete dit;
// Removed 2012-03-08
//	 if (info>1) 
//	    cout<<"The data of the person \""<<p.name()<<"\" has been removed "
//		<<"from system \""<<Systemname<<"\".\n";
	 recalc_data = 1;
	 return 1;
      }
      dit = dit->next;
   }
   return 0;
}


/*
void allelesystem::compute_dataprob() {
  dataitem* dit;
  int i,j;
  for (i=0; i<n_dataalleles; i++) 
  {
     delete[] dataprobmatrixFemale[i];
     delete[] dataprobmatrixMale[i];
  }
  delete[] dataprobmatrixFemale;
  delete[] dataprobmatrixMale;
  delete[] dataprobability;
  delete[] index;
  
  index = new int[n_alleles];
  //The number of alleles in the data matrix that does not 
  //appear in the data (or as a silent allele): 
  int nExtra=0; 

  if (mutationModelFemale < 2 && mutationModelMale < 2)
  {
	//First, find the alleles that appear in the data: 
	int* existsInData = new int[n_alleles]; 
	for (i=0; i<n_alleles; i++) existsInData[i] = 0;
	for (i=0; i<n_alleles; i++) index[i] = 0; 
	dit = data;
	while (dit) {
		existsInData[dit->allele1] = 1;
		existsInData[dit->allele2] = 1;
		dit = dit->next;
	}
	if (hasSilentAllele) 
		  existsInData[silentAllele] = 1;

	//Determine whether an extra allele is needed: 
	double sum = 0; 
	for (i=0; i<n_alleles; i++)
		if (existsInData[i]) sum += probability[i]; 
	if (sum<1) nExtra = 1;  
		
	n_dataalleles = nExtra; 
	//Find dataprobability and index: 
	for (i=0; i<n_alleles; i++) if (existsInData[i]) index[i] = n_dataalleles++;
	dataprobability = new double[n_dataalleles];
	if (nExtra)
	{
		dataprobability[0] = 1;
		for (i=0; i<n_alleles; i++) if (existsInData[i]) 
			dataprobability[0] -= (dataprobability[index[i]] = probability[i]);
	}
	else
	{
		for (i=0; i<n_alleles; i++) if (existsInData[i])
			dataprobability[index[i]] = probability[i];
	}
	delete[] existsInData; 

//	//This may only happen because we have had "approximate" tests above, 
//		//which we do to avoid numerical problems. Thus, this is a numerical 
//	//adjustment: 
//	if (dataprobability[0]<0)
//	{
//		 //Naughty: Renormalization: 
//		for (i=1; i<n_dataalleles; i++)
//		dataprobability[i] /= 1 - dataprobability[0]; 
//		dataprobability[0] = 0; 
//	}

  
  }
	else
	{
		//Determine whether two extra alleles are needed:
		double sum = 0;
		for (i=0; i<n_alleles; i++)
			sum += probability[i];
		if (sum<1) nExtra = 2; 
		n_dataalleles = n_alleles + nExtra; 
		dataprobability = new double[n_dataalleles]; 
		int offset = (sum<1); 
		for (i=0; i<n_alleles; i++)
		{
			index[i] = i+offset; 
			dataprobability[i+offset]=probability[i]; 
		}
		if (sum<1) 
			dataprobability[0]=(dataprobability[n_dataalleles-1]=0.5*(1-sum)); 
	}
 
  //We assume, for both male and female data:
  //0<=mutationrate<1, n_possibilities >= n_dataalleles,
  //and n_possibilities>=2. If n_possibilities==n_dataalleles,
  //then dataprobability[0]==0.

  dataprobmatrixFemale = new double*[n_dataalleles];
  dataprobmatrixMale   = new double*[n_dataalleles];
  for (i=0; i<n_dataalleles; i++) 
  {
     dataprobmatrixFemale[i] = new double[n_dataalleles];
     dataprobmatrixMale[i]   = new double[n_dataalleles];
  }

  // Needed later: 
  // (this "crossum" parameter should be independent of which of the two choices
  // above that have been made)
  double crossum = 0; 
  double totalsum = 0; 
  for (i=0; i<n_alleles; i++) 
  {
	  crossum += probability[i]*(1-probability[i]); 
	  totalsum += probability[i]; 
  }
  crossum += (1-totalsum)*totalsum; //include also the last "extra" allele
  if (crossum == 0) crossum = 1; //avoid numerical problems

  //Given that a mutation happens, there is an equal probability of 
  //ending up in each of the other possible alleles: 
  if (mutationModelFemale==0)
  {
	  if (nExtra==0)
	  {
		  //At this point, we must have n_possibilities==n_alleles==n_dataalleles
		  for (i=0; i<n_dataalleles; i++)
			  for (j=0; j<n_dataalleles; j++)
				  if (i==j)
					  dataprobmatrixFemale[i][j] = 1-mutationrateFemale; 
				  else
					  dataprobmatrixFemale[i][j] = mutationrateFemale/
					  (n_possibilities-1.0); 
	  }
	  else if (nExtra==1)
	  {
		  dataprobmatrixFemale[0][0] = 1.0 - mutationrateFemale*(n_dataalleles-1.0)/
	      (n_possibilities-1.0);
          for (i = 1; i<n_dataalleles; i++) 
	          dataprobmatrixFemale[i][0] = mutationrateFemale*
	          (1.0-(n_dataalleles-2.0)/(n_possibilities-1.0));
          for (i = 0; i<n_dataalleles; i++) 
	          for (j = 1; j<n_dataalleles; j++) 
	              if (i==j) 
	                  dataprobmatrixFemale[i][j] = 1.0 - mutationrateFemale;
	              else
	                  dataprobmatrixFemale[i][j] = mutationrateFemale/
		              (n_possibilities-1.0);
	  }
	  else
	  {
		  for (i=0; i<n_dataalleles; i++)
			  for (j=0; j<n_dataalleles; j++)
				  if ((i==0 && j==0) || (i==n_dataalleles-1 && j==n_dataalleles-1))
					  dataprobmatrixFemale[i][j] = 1-mutationrateFemale +
					  0.5*mutationrateFemale*(n_possibilities-n_dataalleles)/
					  (n_possibilities-1.0); 
				  else if ((i==0 && j==n_dataalleles-1)||(i==n_dataalleles-1 && j==0))
					  dataprobmatrixFemale[i][j] = 
					  0.5*mutationrateFemale*(n_possibilities-n_dataalleles)/
					  (n_possibilities-1.0); 
				  else if (j==0 || j==n_dataalleles-1)
					  dataprobmatrixFemale[i][j] = 
					  0.5*mutationrateFemale*(n_possibilities-n_dataalleles+2.0)/
					  (n_possibilities-1.0); 
				  else if (i==j)
					  dataprobmatrixFemale[i][j] = 
					  1-mutationrateFemale; 
				  else
					  dataprobmatrixFemale[i][j] = 
					  mutationrateFemale/(n_possibilities-1); 
	  }
  }
  //Given that a mutation happens, the probability of ending up in 
  //another allele is proportional to the frequency of that allele: 
  else if (mutationModelFemale==1)
  {
    // Compute alpha: 
	// 
	// Gammelt:
	//
    // double sumsquare = 0; 
    // double sum = 0; 
    // for (i=1; i<n_dataalleles; i++)
    // {
	// sum += dataprobability[i]; 
	// sumsquare += dataprobability[i]*dataprobability[i]; 
    // }
    // sumsquare += (1-sum)*(1-sum)/(n_possibilities-n_dataalleles+1); 
	//
	//double alpha = mutationrateFemale/(1-sumsquare); 
	//
	double alpha = mutationrateFemale/crossum; 

	if (nExtra==0)
	{
		for (i=0; i<n_dataalleles; i++)
			for (j=0; j<n_dataalleles; j++)
				if (i==j)
					dataprobmatrixFemale[i][j] = 1.0 - alpha + alpha*dataprobability[j]; 
				else 
					dataprobmatrixFemale[i][j] = alpha*dataprobability[j]; 
	}
	else if (nExtra==1)
	{
		for (i=0; i<n_dataalleles; i++)
		{
			double tmpsum = 0; 
			for (j=1; j<n_dataalleles; j++)
				if (i==j)
					tmpsum += (dataprobmatrixFemale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
				else
					tmpsum += (dataprobmatrixFemale[i][j] = alpha*dataprobability[j]); 
			dataprobmatrixFemale[i][0] = 1-tmpsum; 
		}
	}
	else
	{
		for (i=0; i<n_dataalleles; i++)
		{
			double tmpsum = 0; 
			for (j=1; j<n_dataalleles-1; j++)
				if (i==j)
					tmpsum += (dataprobmatrixFemale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
				else
					tmpsum += (dataprobmatrixFemale[i][j] = alpha*dataprobability[j]); 
			dataprobmatrixFemale[i][0] = dataprobmatrixFemale[i][n_dataalleles-1] = 0.5*(1-tmpsum); 
		}
	}
  }
  //Given that a mutation happens, the ratio of probabilities of ending up in 
  //two different alleles is a constant to the power of the difference of
  //their distances to the start allele:
  else if (mutationModelFemale==2)
  {
	  for (i=0; i<n_dataalleles; i++)
      {
    	  double k = mutationrateFemale*(1-mutationRangeFemale)/
	          mutationRangeFemale/(2-mypow(mutationRangeFemale, i)
	          -mypow(mutationRangeFemale, n_dataalleles-i-1)); 
		  for (j=0; j<n_dataalleles; j++)
	          if (i==j)
	              dataprobmatrixFemale[i][j] = 1-mutationrateFemale;  
	          else if (j<i)
	              dataprobmatrixFemale[i][j] = k*mypow(mutationRangeFemale, i-j); 
	          else
	              dataprobmatrixFemale[i][j] = k*mypow(mutationRangeFemale, j-i); 
	  }


//     //For simplicity, we first set up a full transition matrix for the
//     //full set of alleles, and then derive the dataprobmatrix from it. 
//	
//	//Find out how many alleles to use: 
//	  double probsum = 0; 
//		for (i=0; i<n_alleles; i++) probsum += probability[i]; 
//	  int nAllelesInRow = n_alleles + 2*(probsum<1) - hasSilentAllele;
//	  
//	  double** tmpMatrix = new double*[nAllelesInRow]; 
//	  for (i=0; i<nAllelesInRow; i++)
//	  {
//		  tmpMatrix[i] = new double[nAllelesInRow]; 
//		  double ki = mutationRateFemale*(1-mutationRangeFemale)/
//			  mutationRangeFemale/(2.0-mypow(mutationRangeFemale, i-1)
//			  -mypow(mutationRangeFemale, nAllelesInRow)); 
//		  if (hasSilentAllele) ki *= (1-1/nAllelesInRow); 
//		  for (j=0; j<nAllelesInRow)
//			  if (i<j)
//				tmpMatrix[i][j] = ki*mypow(mutationRangeFemale, j-i)
//			  else if (i>j)
//			    tmpMatrix[i][j] = ki*mypow(mutationRangeFemale, i-j)
//			  else
//			    tmpMatrix[i][j] = 1-mutationRateFemale; 
//	  }


//     //For now: 
//     //For even more simplicity (and accuracy) we just use the full matrix
//     delete[] dataprobability; 
//     for (i=0; i<n_dataalleles; i++) 
//	 delete[] dataprobmatrixFemale[i]; 
//     delete[] dataprobmatrixFemale; 
//
//     double sum = 0; 
//     for (i=0; i<n_alleles; i++)
//	 sum += probability[i]; 
//
//     if (sum>=1)
//     {
//	n_dataalleles = n_alleles; 
//	dataprobability = new double[n_dataalleles]; 
//	for (i=0; i<n_alleles; i++)
//	{
//	   index[i] = i; 
//	   dataprobability[i] = probability[i]; 
//	}
//
//}
//     else
//     {
//	n_dataalleles = n_alleles+2; 
//	dataprobability = new double[n_dataalleles]; 
//	dataprobability[n_dataalleles-1] = 
//	   dataprobability[0] = 0.5*(1-sum); 
//	for (i=0; i<n_alleles; i++)
//	{
//	   index[i] = i+1; 
//	   dataprobability[i+1] = probability[i]; 
//	}
//     }
//
//     dataprobmatrixFemale  = new double*[n_dataalleles]; 
//     for (i=0; i<n_dataalleles; i++)
//     {
//	dataprobmatrixFemale[i]  = new double[n_dataalleles]; 
//	double k = mutationrateFemale*(1-mutationRangeFemale)/
//	   mutationRangeFemale/
//	   (2-mypow(mutationRangeFemale, i-1)
//	    -mypow(mutationRangeFemale, n_dataalleles-i)); 
//	for (j=0; j<n_dataalleles; j++)
//	   if (i==j)
//	      dataprobmatrixFemale[i][j] = 1-mutationrateFemale; 
//	   else if (j<i)
//	      dataprobmatrixFemale[i][j] = k*mypow(mutationRangeFemale, i-j); 
//	   else
//	      dataprobmatrixFemale[i][j] = k*mypow(mutationRangeFemale, j-i); 
//     }
  }
  else //mutationModel==3
  {
	  double constant = mutationrateFemale*(1-mutationRangeFemale)*(1-mutationRangeFemale)/
		  2.0/mutationRangeFemale/(n_dataalleles-mutationRangeFemale*n_dataalleles-1+
		  mypow(mutationRangeFemale, n_dataalleles)); 
	  double sum = 0; 
	  for (i=0; i<n_dataalleles; i++)
      {
		  for (j=0; j<n_dataalleles; j++)
			  if (i<j) {
	              dataprobmatrixFemale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeFemale, j-i); 
				  sum += dataprobmatrixFemale[i][j]; 
			  }
			  else if (j<i) {
	              dataprobmatrixFemale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeFemale, i-j);
				  sum += dataprobmatrixFemale[i][j]; 
			  }
		  dataprobmatrixFemale[i][i] = 1-sum; 
	  }
  }


  //Given that a mutation happens, there is an equal probability of 
  //ending up in each of the other possible alleles: 
  if (mutationModelMale==0)
  {
	  if (nExtra==0)
	  {
		  //At this point, we must have n_possibilities==n_alleles==n_dataalleles
		  for (i=0; i<n_dataalleles; i++)
			  for (j=0; j<n_dataalleles; j++)
				  if (i==j)
					  dataprobmatrixMale[i][j] = 1-mutationrateMale; 
				  else
					  dataprobmatrixMale[i][j] = mutationrateMale/
					  (n_possibilities-1.0); 
	  }
	  else if (nExtra==1)
	  {
		  dataprobmatrixMale[0][0] = 1.0 - mutationrateMale*(n_dataalleles-1.0)/
	      (n_possibilities-1.0);
          for (i = 1; i<n_dataalleles; i++) 
	          dataprobmatrixMale[i][0] = mutationrateMale*
	          (1.0-(n_dataalleles-2.0)/(n_possibilities-1.0));
          for (i = 0; i<n_dataalleles; i++) 
	          for (j = 1; j<n_dataalleles; j++) 
	              if (i==j) 
	                  dataprobmatrixMale[i][j] = 1.0 - mutationrateMale;
	              else
	                  dataprobmatrixMale[i][j] = mutationrateMale/
		              (n_possibilities-1.0);
	  }
	  else
	  {
		  for (i=0; i<n_dataalleles; i++)
			  for (j=0; j<n_dataalleles; j++)
				  if ((i==0 && j==0) || (i==n_dataalleles-1 && j==n_dataalleles-1))
					  dataprobmatrixMale[i][j] = 1-mutationrateMale +
					  0.5*mutationrateMale*(n_possibilities-n_dataalleles)/
					  (n_possibilities-1.0); 
				  else if ((i==0 && j==n_dataalleles-1)||(i==n_dataalleles-1 && j==0))
					  dataprobmatrixMale[i][j] = 
					  0.5*mutationrateMale*(n_possibilities-n_dataalleles)/
					  (n_possibilities-1.0); 
				  else if (j==0 || j==n_dataalleles-1)
					  dataprobmatrixMale[i][j] = 
					  0.5*mutationrateMale*(n_possibilities-n_dataalleles+2.0)/
					  (n_possibilities-1.0); 
				  else if (i==j)
					  dataprobmatrixMale[i][j] = 
					  1-mutationrateMale; 
				  else
					  dataprobmatrixMale[i][j] = 
					  mutationrateMale/(n_possibilities-1); 
	  }
  }
  //Given that a mutation happens, the probability of ending up in 
  //another allele is proportional to the frequency of that allele: 
  else if (mutationModelMale==1)
  {
	double alpha = mutationrateMale/crossum; 
	if (nExtra==0)
	{
		for (i=0; i<n_dataalleles; i++)
			for (j=0; j<n_dataalleles; j++)
				if (i==j)
					dataprobmatrixMale[i][j] = 1.0 - alpha + alpha*dataprobability[j]; 
				else 
					dataprobmatrixMale[i][j] = alpha*dataprobability[j]; 
	}
	else if (nExtra==1)
	{
		for (i=0; i<n_dataalleles; i++)
		{
			double tmpsum = 0; 
			for (j=1; j<n_dataalleles; j++)
				if (i==j)
					tmpsum += (dataprobmatrixMale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
				else
					tmpsum += (dataprobmatrixMale[i][j] = alpha*dataprobability[j]); 
			dataprobmatrixMale[i][0] = 1-tmpsum; 
		}
	}
	else
	{
		for (i=0; i<n_dataalleles; i++)
		{
			double tmpsum = 0; 
			for (j=1; j<n_dataalleles-1; j++)
				if (i==j)
					tmpsum += (dataprobmatrixMale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
				else
					tmpsum += (dataprobmatrixMale[i][j] = alpha*dataprobability[j]); 
			dataprobmatrixMale[i][0] = dataprobmatrixMale[i][n_dataalleles-1] = 0.5*(1-tmpsum); 
		}
	}
  }
  //Given that a mutation happens, the ratio of probabilities of ending up in 
  //two different alleles is a constant to the power of the difference of
  //their distances to the start allele:
  else if (mutationModelMale==2)
  {
	  for (i=0; i<n_dataalleles; i++)
      {
    	  double k = mutationrateMale*(1-mutationRangeMale)/
	          mutationRangeMale/(2-mypow(mutationRangeMale, i)
	          -mypow(mutationRangeMale, n_dataalleles-i-1)); 
	      for (j=0; j<n_dataalleles; j++)
	          if (i==j)
	              dataprobmatrixMale[i][j] = 1-mutationrateMale;  
	          else if (j<i)
	              dataprobmatrixMale[i][j] = k*mypow(mutationRangeMale, i-j); 
	          else
	              dataprobmatrixMale[i][j] = k*mypow(mutationRangeMale, j-i); 
	  }
  }
  else //mutationModel==3
  {
	  double constant = mutationrateMale*(1-mutationRangeMale)*(1-mutationRangeMale)/
		  2.0/mutationRangeMale/(n_dataalleles-mutationRangeMale*n_dataalleles-1+
		  mypow(mutationRangeMale, n_dataalleles)); 
	  double sum = 0; 
	  for (i=0; i<n_dataalleles; i++)
      {
		  for (j=0; j<n_dataalleles; j++)
			  if (i<j) {
	              dataprobmatrixMale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeMale, j-i); 
				  sum += dataprobmatrixMale[i][j]; 
			  }
			  else if (j<i) {
	              dataprobmatrixMale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeMale, i-j);
				  sum += dataprobmatrixMale[i][j]; 
			  }
		  dataprobmatrixMale[i][i] = 1-sum; 
	  }
  }
*/


/*
  //Given that a mutation happens, there is an equal probability of 
  //ending up in each of the other possible alleles: 
  if (mutationModelMale==0)
  {
     dataprobmatrixMale[0][0] = 1.0 - mutationrateMale*(n_dataalleles-1.0)/
	(n_possibilities-1.0);
     for (i = 1; i<n_dataalleles; i++) 
	dataprobmatrixMale[i][0] = mutationrateMale*
	   (1.0-(n_dataalleles-2.0)/(n_possibilities-1.0));
     for (i = 0; i<n_dataalleles; i++) 
	for (j = 1; j<n_dataalleles; j++) 
	   if (i==j) 
	      dataprobmatrixMale[i][j] = 1.0 - mutationrateMale;
	   else
	      dataprobmatrixMale[i][j] = mutationrateMale/
		 (n_possibilities-1.0);
  }
  //Given that a mutation happens, the probability of ending up in 
  //another allele is proportional to the frequency of that allele: 
  else if (mutationModelMale==1)
  {
     //Compute alpha: 
     //
	 // GAMMELT:
	 // double sumsquare = 0; 
     // double sum = 0; 
     // for (i=1; i<n_dataalleles; i++)
     // {
	 // sum += dataprobability[i]; 
	 // sumsquare += dataprobability[i]*dataprobability[i]; 
     // }
     // sumsquare += (1-sum)*(1-sum)/(n_possibilities-n_dataalleles+1); 
     // double alpha = mutationrateMale/(1-sumsquare); 
	 //
	 double alpha = mutationrateMale/crossum; 

     for (i=0; i<n_dataalleles; i++)
	for (j=0; j<n_dataalleles; j++)
	   if (i==j)
	      dataprobmatrixMale[i][j] = 1.0 - alpha + alpha*dataprobability[j]; 
	   else
	      dataprobmatrixMale[i][j] = alpha*dataprobability[j]; 
  }
  //Given that a mutation happens, the ratio of probabilities of ending up in 
  //two different alleles is a constant to the power of the difference of
  //their distances to the start allele:
  else //mutationModel==2
  {
     //For simplicity, we first set up a full transition matrix for the
     //full set of alleles, and then derive the dataprobmatrix from it. 
     

     //For now: 
     //For even more simplicity (and accuracy) we just use the full matrix
     delete[] dataprobability; 
     for (i=0; i<n_dataalleles; i++) 
	delete[] dataprobmatrixMale[i]; 
     delete[] dataprobmatrixMale; 

     double sum = 0; 
     for (i=0; i<n_alleles; i++)
	sum += probability[i]; 

     if (sum>=1)
     {
	n_dataalleles = n_alleles; 
	dataprobability = new double[n_dataalleles]; 
	for (i=0; i<n_alleles; i++)
	{
	   index[i] = i; 
	   dataprobability[i] = probability[i]; 
	}

     }
     else
     {
	n_dataalleles = n_alleles+2; 
	dataprobability = new double[n_dataalleles]; 
	dataprobability[n_dataalleles-1] = 
	   dataprobability[0] = 0.5*(1-sum); 
	for (i=0; i<n_alleles; i++)
	{
	   index[i] = i+1; 
	   dataprobability[i+1] = probability[i]; 
	}
     }

     dataprobmatrixMale  = new double*[n_dataalleles]; 
     for (i=0; i<n_dataalleles; i++)
     {
	dataprobmatrixMale[i]  = new double[n_dataalleles]; 
	double k = mutationrateMale*(1-mutationRangeMale)/
	   mutationRangeMale/
	   (2-mypow(mutationRangeMale, i-1)
	    -mypow(mutationRangeMale, n_dataalleles-i)); 
	for (j=0; j<n_dataalleles; j++)
	   if (i==j)
	      dataprobmatrixMale[i][j] = 1-mutationrateMale; 
	   else if (j<i)
	      dataprobmatrixMale[i][j] = k*mypow(mutationRangeMale, i-j); 
	   else
	      dataprobmatrixMale[i][j] = k*mypow(mutationRangeMale, j-i); 
     }
  }
*/

/*
  //Midlertidig utskrift av transisjonsmatrise: 
  //printMutMatrix(); 

  recalc_data = 0;
}

*/



//New version, after deciding that all systems should have sum of frequencies equal to 1:
void allelesystem::compute_dataprob() {
  dataitem* dit;
  int i,j;
  for (i=0; i<n_dataalleles; i++) 
  {
     delete[] dataprobmatrixFemale[i];
     delete[] dataprobmatrixMale[i];
  }
  delete[] dataprobmatrixFemale;
  delete[] dataprobmatrixMale;
  delete[] dataprobability;
  delete[] index;
  

  //THIS IS THE NEW VERSION; WHERE WE ALWAYS HAVE
  n_possibilities = n_alleles; 


  index = new int[n_alleles];
  //The number of alleles in the data matrix that does not 
  //appear in the data (or as a silent allele): 
  int nExtra=0; 

///////////////////////////////////////////////////
// LAST ADDITION (2004-08-30): Change to make 
// the program run faster when the mutation rate is zero, 
// or the system has no data: 
///////////////////////////////////////////////////

  int femaleMutModelUsed = mutationModelFemale; 
  int maleMutModelUsed   = mutationModelMale; 
  if (mutationrateFemale == 0 && mutationrateMale == 0) 
  {
	  femaleMutModelUsed = 0; 
	  maleMutModelUsed   = 0; 
  }
  if (data == 0) 
  {
	  femaleMutModelUsed = 0; 
	  maleMutModelUsed  = 0; 
  }



///////////////////////////////////////////////////

  if (femaleMutModelUsed < 2 && maleMutModelUsed < 2)
  {
	//First, find the alleles that appear in the data: 
	int* existsInData = new int[n_alleles]; 
	for (i=0; i<n_alleles; i++) existsInData[i] = 0;
	for (i=0; i<n_alleles; i++) index[i] = 0; 
	dit = data;
	while (dit) {
		existsInData[dit->allele1] = 1;
		existsInData[dit->allele2] = 1;
		dit = dit->next;
	}
	if (hasSilentAllele) 
		  existsInData[silentAllele] = 1;

	//Determine whether an extra allele is needed: 
//	double sum = 0; 
//	for (i=0; i<n_alleles; i++)
//		if (existsInData[i]) sum += probability[i]; 
//	if (sum<1) nExtra = 1;  
	//ALTERNATIVE TO ABOVE: 
	for (i=0; i<n_alleles; i++)
		if (!existsInData[i]) nExtra = 1; 
		
	n_dataalleles = nExtra; 
	//Find dataprobability and index: 
	for (i=0; i<n_alleles; i++) if (existsInData[i]) index[i] = n_dataalleles++;
	dataprobability = new double[n_dataalleles];
	if (nExtra)
	{
		dataprobability[0] = 1;
		for (i=0; i<n_alleles; i++) if (existsInData[i]) 
			dataprobability[0] -= (dataprobability[index[i]] = probability[i]);
	}
	else
	{
		for (i=0; i<n_alleles; i++) if (existsInData[i])
			dataprobability[index[i]] = probability[i];
	}
	delete[] existsInData; 
  
  }
	else
	{
		//Determine whether two extra alleles are needed:
//		double sum = 0;
//		for (i=0; i<n_alleles; i++)
//			sum += probability[i];
//		if (sum<1) nExtra = 2; 
		n_dataalleles = n_alleles + nExtra; 
		dataprobability = new double[n_dataalleles]; 
		int offset = 0; 
//		int offset = (sum<1); 
		for (i=0; i<n_alleles; i++)
		{
			index[i] = i+offset; 
			dataprobability[i+offset]=probability[i]; 
		}
//		if (sum<1) 
//			dataprobability[0]=(dataprobability[n_dataalleles-1]=0.5*(1-sum)); 
	}
 
  //We assume, for both male and female data:
  //0<=mutationrate<1, n_possibilities >= n_dataalleles,
  //and n_possibilities>=2. If n_possibilities==n_dataalleles,
  //then dataprobability[0]==0.

  dataprobmatrixFemale = new double*[n_dataalleles];
  dataprobmatrixMale   = new double*[n_dataalleles];
  for (i=0; i<n_dataalleles; i++) 
  {
     dataprobmatrixFemale[i] = new double[n_dataalleles];
     dataprobmatrixMale[i]   = new double[n_dataalleles];
  }

  // Needed later: 
  // (this "crossum" parameter should be independent of which of the two choices
  // above that have been made)
  double crossum = 0; 
  double totalsum = 0; 
  for (i=0; i<n_alleles; i++) 
  {
	  crossum += probability[i]*(1-probability[i]); 
	  totalsum += probability[i]; 
  }
//  crossum += (1-totalsum)*totalsum; //include also the last "extra" allele
  if (crossum == 0) crossum = 1; //avoid numerical problems

  //Given that a mutation happens, there is an equal probability of 
  //ending up in each of the other possible alleles: 
  if (femaleMutModelUsed==0)
  {
	  if (nExtra==0)
	  {
		  //At this point, we must have n_possibilities==n_alleles==n_dataalleles
		  for (i=0; i<n_dataalleles; i++)
			  for (j=0; j<n_dataalleles; j++)
				  if (i==j)
					  dataprobmatrixFemale[i][j] = 1-mutationrateFemale; 
				  else
					  dataprobmatrixFemale[i][j] = mutationrateFemale/
					  (n_possibilities-1.0); 
	  }
	  else //if (nExtra==1)
	  {
		  dataprobmatrixFemale[0][0] = 1.0 - mutationrateFemale*(n_dataalleles-1.0)/
	      (n_possibilities-1.0);
          for (i = 1; i<n_dataalleles; i++) 
	          dataprobmatrixFemale[i][0] = mutationrateFemale*
	          (1.0-(n_dataalleles-2.0)/(n_possibilities-1.0));
          for (i = 0; i<n_dataalleles; i++) 
	          for (j = 1; j<n_dataalleles; j++) 
	              if (i==j) 
	                  dataprobmatrixFemale[i][j] = 1.0 - mutationrateFemale;
	              else
	                  dataprobmatrixFemale[i][j] = mutationrateFemale/
		              (n_possibilities-1.0);
	  }
//	  else
//	  {
//		  for (i=0; i<n_dataalleles; i++)
//			  for (j=0; j<n_dataalleles; j++)
//				  if ((i==0 && j==0) || (i==n_dataalleles-1 && j==n_dataalleles-1))
//					  dataprobmatrixFemale[i][j] = 1-mutationrateFemale +
//					  0.5*mutationrateFemale*(n_possibilities-n_dataalleles)/
//					  (n_possibilities-1.0); 
//				  else if ((i==0 && j==n_dataalleles-1)||(i==n_dataalleles-1 && j==0))
//					  dataprobmatrixFemale[i][j] = 
//					  0.5*mutationrateFemale*(n_possibilities-n_dataalleles)/
//					  (n_possibilities-1.0); 
//				  else if (j==0 || j==n_dataalleles-1)
//					  dataprobmatrixFemale[i][j] = 
//					  0.5*mutationrateFemale*(n_possibilities-n_dataalleles+2.0)/
//					  (n_possibilities-1.0); 
//				  else if (i==j)
//					  dataprobmatrixFemale[i][j] = 
//					  1-mutationrateFemale; 
//				  else
//					  dataprobmatrixFemale[i][j] = 
//					  mutationrateFemale/(n_possibilities-1); 
//	  }
  }
  //Given that a mutation happens, the probability of ending up in 
  //another allele is proportional to the frequency of that allele: 
  else if (femaleMutModelUsed==1)
  {
	double alpha = mutationrateFemale/crossum; 

	if (nExtra==0)
	{
		for (i=0; i<n_dataalleles; i++)
			for (j=0; j<n_dataalleles; j++)
				if (i==j)
					dataprobmatrixFemale[i][j] = 1.0 - alpha + alpha*dataprobability[j]; 
				else 
					dataprobmatrixFemale[i][j] = alpha*dataprobability[j]; 
	}
	else //if (nExtra==1)
	{
		for (i=0; i<n_dataalleles; i++)
		{
			double tmpsum = 0; 
			for (j=1; j<n_dataalleles; j++)
				if (i==j)
					tmpsum += (dataprobmatrixFemale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
				else
					tmpsum += (dataprobmatrixFemale[i][j] = alpha*dataprobability[j]); 
			dataprobmatrixFemale[i][0] = 1-tmpsum; 
		}
	}
//	else
//	{
//		for (i=0; i<n_dataalleles; i++)
//		{
//			double tmpsum = 0; 
//			for (j=1; j<n_dataalleles-1; j++)
//				if (i==j)
//					tmpsum += (dataprobmatrixFemale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
//				else
//					tmpsum += (dataprobmatrixFemale[i][j] = alpha*dataprobability[j]); 
//			dataprobmatrixFemale[i][0] = dataprobmatrixFemale[i][n_dataalleles-1] = 0.5*(1-tmpsum); 
//		}
//	}
  }
  //Given that a mutation happens, the ratio of probabilities of ending up in 
  //two different alleles is a constant to the power of the difference of
  //their distances to the start allele:
  else if (femaleMutModelUsed==2)
  {
	  for (i=0; i<n_dataalleles; i++)
      {
    	  double k = mutationrateFemale*(1-mutationRangeFemale)/
	          mutationRangeFemale/(2-mypow(mutationRangeFemale, i)
	          -mypow(mutationRangeFemale, n_dataalleles-i-1)); 
		  for (j=0; j<n_dataalleles; j++)
	          if (i==j)
	              dataprobmatrixFemale[i][j] = 1-mutationrateFemale;  
	          else if (j<i)
	              dataprobmatrixFemale[i][j] = k*mypow(mutationRangeFemale, i-j); 
	          else
	              dataprobmatrixFemale[i][j] = k*mypow(mutationRangeFemale, j-i); 
	  }

  }
  else //mutationModel==3
  {
	  double constant = mutationrateFemale*(1-mutationRangeFemale)*(1-mutationRangeFemale)/
		  2.0/mutationRangeFemale/(n_dataalleles-mutationRangeFemale*n_dataalleles-1+
		  mypow(mutationRangeFemale, n_dataalleles)); 
	  for (i=0; i<n_dataalleles; i++)
      {
		  double sum = 0; 
		  for (j=0; j<n_dataalleles; j++)
			  if (i<j) {
	              dataprobmatrixFemale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeFemale, j-i); 
				  sum += dataprobmatrixFemale[i][j]; 
			  }
			  else if (j<i) {
	              dataprobmatrixFemale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeFemale, i-j);
				  sum += dataprobmatrixFemale[i][j]; 
			  }
		  dataprobmatrixFemale[i][i] = 1-sum; 
	  }
  }


  //Given that a mutation happens, there is an equal probability of 
  //ending up in each of the other possible alleles: 
  if (maleMutModelUsed==0)
  {
	  if (nExtra==0)
	  {
		  //At this point, we must have n_possibilities==n_alleles==n_dataalleles
		  for (i=0; i<n_dataalleles; i++)
			  for (j=0; j<n_dataalleles; j++)
				  if (i==j)
					  dataprobmatrixMale[i][j] = 1-mutationrateMale; 
				  else
					  dataprobmatrixMale[i][j] = mutationrateMale/
					  (n_possibilities-1.0); 
	  }
	  else //if (nExtra==1)
	  {
		  dataprobmatrixMale[0][0] = 1.0 - mutationrateMale*(n_dataalleles-1.0)/
	      (n_possibilities-1.0);
          for (i = 1; i<n_dataalleles; i++) 
	          dataprobmatrixMale[i][0] = mutationrateMale*
	          (1.0-(n_dataalleles-2.0)/(n_possibilities-1.0));
          for (i = 0; i<n_dataalleles; i++) 
	          for (j = 1; j<n_dataalleles; j++) 
	              if (i==j) 
	                  dataprobmatrixMale[i][j] = 1.0 - mutationrateMale;
	              else
	                  dataprobmatrixMale[i][j] = mutationrateMale/
		              (n_possibilities-1.0);
	  }
//	  else
//	  {
//		  for (i=0; i<n_dataalleles; i++)
//			  for (j=0; j<n_dataalleles; j++)
//				  if ((i==0 && j==0) || (i==n_dataalleles-1 && j==n_dataalleles-1))
//					  dataprobmatrixMale[i][j] = 1-mutationrateMale +
//					  0.5*mutationrateMale*(n_possibilities-n_dataalleles)/
//					  (n_possibilities-1.0); 
//				  else if ((i==0 && j==n_dataalleles-1)||(i==n_dataalleles-1 && j==0))
//					  dataprobmatrixMale[i][j] = 
//					  0.5*mutationrateMale*(n_possibilities-n_dataalleles)/
//					  (n_possibilities-1.0); 
//				  else if (j==0 || j==n_dataalleles-1)
//					  dataprobmatrixMale[i][j] = 
//					  0.5*mutationrateMale*(n_possibilities-n_dataalleles+2.0)/
//					  (n_possibilities-1.0); 
//				  else if (i==j)
//					  dataprobmatrixMale[i][j] = 
//					  1-mutationrateMale; 
//				  else
//					  dataprobmatrixMale[i][j] = 
//					  mutationrateMale/(n_possibilities-1); 
//	  }
  }
  //Given that a mutation happens, the probability of ending up in 
  //another allele is proportional to the frequency of that allele: 
  else if (maleMutModelUsed==1)
  {
	double alpha = mutationrateMale/crossum; 
	if (nExtra==0)
	{
		for (i=0; i<n_dataalleles; i++)
			for (j=0; j<n_dataalleles; j++)
				if (i==j)
					dataprobmatrixMale[i][j] = 1.0 - alpha + alpha*dataprobability[j]; 
				else 
					dataprobmatrixMale[i][j] = alpha*dataprobability[j]; 
	}
	else //if (nExtra==1)
	{
		for (i=0; i<n_dataalleles; i++)
		{
			double tmpsum = 0; 
			for (j=1; j<n_dataalleles; j++)
				if (i==j)
					tmpsum += (dataprobmatrixMale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
				else
					tmpsum += (dataprobmatrixMale[i][j] = alpha*dataprobability[j]); 
			dataprobmatrixMale[i][0] = 1-tmpsum; 
		}
	}
//	else
//	{
//		for (i=0; i<n_dataalleles; i++)
//		{
//			double tmpsum = 0; 
//			for (j=1; j<n_dataalleles-1; j++)
//				if (i==j)
//					tmpsum += (dataprobmatrixMale[i][j] = 1.0 - alpha + alpha*dataprobability[j]); 
//				else
//					tmpsum += (dataprobmatrixMale[i][j] = alpha*dataprobability[j]); 
//			dataprobmatrixMale[i][0] = dataprobmatrixMale[i][n_dataalleles-1] = 0.5*(1-tmpsum); 
//		}
//	}
  }
  //Given that a mutation happens, the ratio of probabilities of ending up in 
  //two different alleles is a constant to the power of the difference of
  //their distances to the start allele:
  else if (maleMutModelUsed==2)
  {
	  for (i=0; i<n_dataalleles; i++)
      {
    	  double k = mutationrateMale*(1-mutationRangeMale)/
	          mutationRangeMale/(2-mypow(mutationRangeMale, i)
	          -mypow(mutationRangeMale, n_dataalleles-i-1)); 
	      for (j=0; j<n_dataalleles; j++)
	          if (i==j)
	              dataprobmatrixMale[i][j] = 1-mutationrateMale;  
	          else if (j<i)
	              dataprobmatrixMale[i][j] = k*mypow(mutationRangeMale, i-j); 
	          else
	              dataprobmatrixMale[i][j] = k*mypow(mutationRangeMale, j-i); 
	  }
  }
  else //mutationModel==3
  {
	  double constant = mutationrateMale*(1-mutationRangeMale)*(1-mutationRangeMale)/
		  2.0/mutationRangeMale/(n_dataalleles-mutationRangeMale*n_dataalleles-1+
		  mypow(mutationRangeMale, n_dataalleles)); 
	  for (i=0; i<n_dataalleles; i++)
      {
		  double sum = 0; 
		  for (j=0; j<n_dataalleles; j++)
			  if (i<j) {
	              dataprobmatrixMale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeMale, j-i); 
				  sum += dataprobmatrixMale[i][j]; 
			  }
			  else if (j<i) {
	              dataprobmatrixMale[i][j] = constant/dataprobability[i]*
				  mypow(mutationRangeMale, i-j);
				  sum += dataprobmatrixMale[i][j]; 
			  }
		  dataprobmatrixMale[i][i] = 1-sum; 
	  }
  }


  //Midlertidig utskrift av transisjonsmatrise: 
  //printMutMatrix(); 

  recalc_data = 0;
}





void allelesystem::printMutMatrix()
{
   int i; 
   int* invindex = new int[n_dataalleles]; 
   for (i=0; i<n_dataalleles; i++)
      invindex[i] = -1; 
   for (i=0; i<n_alleles; i++)
      if (invindex[index[i]]>=0)
	 invindex[index[i]] = 0; 
      else
	 invindex[index[i]] = i+1; 
   
   ofstream of1("MutModelFemale.txt"); 
   ofstream of2("MutModelMale.txt"); 
   if (of1.good() && of2.good())
   {
      of1<<"     "; 
      of2<<"     "; 
      for (i=0; i<n_dataalleles; i++)
      {
	 of1.width(12); 
	 of2.width(12); 
	 of1<<invindex[i]; 
	 of2<<invindex[i]; 
      }
      of1<<'\n'; 
      of2<<'\n'; 
      for (i=0; i<n_dataalleles; i++)
      {
	 of1.width(3); 
	 of2.width(3); 
	 of1<<invindex[i]<<"  "; 
	 of2<<invindex[i]<<"  "; 
	 for (int j=0; j<n_dataalleles; j++)
	 {
	    of1.width(12); 
	    of2.width(12); 
	    of1<<dataprobmatrixFemale[i][j]; 
	    of2<<dataprobmatrixMale[i][j]; 
	 }
	 of1<<'\n'; 
	 of2<<'\n'; 
      }
   }
   delete[] invindex; 
}

void allelesystem::execute(family& fam, 
			   int info, 
			   int& error) 
{
   if (recalc_data) compute_dataprob();

   systemdata sd(Systemname, n_dataalleles,
		 dataprobability, dataprobmatrixFemale, 
		 dataprobmatrixMale, kinship,
		 hasSilentAllele, 
		 (hasSilentAllele)? index[silentAllele] : 0);
   dataitem* dit = data;
   while (dit) {
      if (fam.add_data(sd, dit->p, index[dit->allele1], 
		       index[dit->allele2], 
		       info, error)) 
      {
	 //The odds collapse is incompatible with the data.
	 result = 0;
	 fam.remove_data();
	 return;
      }
      dit = dit->next;
   }
   result = fam.execute(sd, error);
   if (error) 
   {
// Removed 2012-03-08
//      if (info>0)
//	 cout<<"ERROR: Too many people in some cutsets for allele system "
//	     <<Systemname<<".\n";
   } 
   else 
   {
// Removed 2012-03-08
//      if (info>1) 
//	 cout<<"Finished computations for system "<<Systemname<<".\n";
   }
   fam.remove_data();
}

void allelesystem::write_freq(ostream& out) {
    out<<separator<<"ALLELE SYSTEM "<<Systemname<<'\n'<<separator;
    out<<"\nMutation probability: "<<0.5*(mutationrateFemale+mutationrateMale);
//    if (mutationrateFemale+mutationrateMale>0)
//	out<<", number of possible alleles: "<<n_possibilities<<"\n\n";
//    else out<<"\n\n";         
	out<<"\n\n";

    if (n_alleles) {
	out<<"General population frequencies of alleles:\n";
	table tab("allele","frequency");
	for (int i=0; i<n_alleles; i++) {
	    tab.put(name[i]);
	    tab.endcolumn();
	    tab.put(probability[i]);
	    tab.endcolumn();
	}
	tab.printout(out);
    } else
	out<<"No alleles registered.\n";
}

void allelesystem::write(ostream& out, int old_results_valid,
			 oddsobject* odds) {
    write_freq(out);
    if (data) {
	out<<"\nObserved alleles in this system:\n";
	table tab("person","observed alleles");
	dataitem* dit = data;
	do {
	    tab.put(dit->p->name());
	    tab.endcolumn();
	    tab.put(name[dit->allele1]);
	    tab.put(name[dit->allele2]);
	    tab.endcolumn();
	} while ((dit = dit->next));
	tab.printout(out);
    } else
	out<<"\nNo observations of alleles registered.\n";
    if (old_results_valid) {
	if (odds) 
	    out<<"\nThe odds that "<<odds->pers1->name()
	       <<" = "<<odds->pers2->name()<<": "<<result<<"\n";
	else
	    out<<"\nThe probability of the data given the family "
	       <<"structure: "<<result<<'\n';
    }
}


