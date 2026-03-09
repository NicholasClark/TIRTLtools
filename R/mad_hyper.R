#madhyper
estimate_pair_prob<-function(wi,wj,w_ij,w_tot,cpw,alpha,prior=1)
{
  fr_match=estimate_match_frequencies(wi,wj,w_ij,w_tot,cpw,alpha)
  fr_nonmatch=estimate_nonmatch_frequencies(wi,wj,w_ij,w_tot,cpw,alpha)
  match_probs=estimate_probs_vec(wi,wj,w_ij,w_tot,cpw,fr_match)
  nonmatch_probs=estimate_probs_vec(wi,wj,w_ij,w_tot,cpw,fr_nonmatch)
  prior*match_probs/nonmatch_probs
}

well_freq<-function(f,c){
  1-(1-f)**c
}

multinom_coeff<-function(vec)
{
  factorial(sum(vec))/prod(factorial(vec)) 
}

estimate_probs_vec<-function(wi,wj,w_ij,w_tot,cpw,freqvec)
{ 
  p_total=1
  for (i in 1:length(wi)){
  wo=w_tot[i]-(wi[i]+wj[i]+w_ij[i])
  fs<-well_freq(freqvec,c = cpw[i])
  p_o = (1-fs[1])*(1-fs[2])*(1-fs[3])
  p_i = fs[1]*(1-fs[2])*(1-fs[3])
  p_j = fs[2]*(1-fs[1])*(1-fs[3])
  p_ij = 1-p_o-p_i-p_j
  p_total=p_total*dmultinom(c(wi[i],wj[i],w_ij[i],wo),prob = c(p_i, p_j, p_ij, p_o))
  }
  p_total
}

estimate_nonmatch_frequencies<-function(wi,wj,w_ij,w_tot,cpw,alpha){
  freqs_i=find_freq(wi+w_ij,w_tot,cpw,alpha)
  freqs_j=find_freq(wj+w_ij,w_tot,cpw,alpha)
  freqs_ij=0
  c(freqs_i,freqs_j,freqs_ij)
}

estimate_match_frequencies<-function(wi,wj,w_ij,w_tot,cpw,alpha){
  w_o=w_tot-(wi+wj+w_ij)
  freqs_i=find_freq(wi,wi+w_o,cpw,alpha) # why it is wi +wo? 
  freqs_j=find_freq(wj,wj+w_o,cpw,alpha) # why it is wi +wo? 
  wij_clonal=wij_adjustment(w_ij,w_tot,cpw,freqs_i,freqs_j)
  freqs_ij<-find_freq(w = wij_clonal,w_tot = (w_tot-w_ij+wij_clonal),cpw = cpw,a = alpha)
  c(freqs_i,freqs_j,freqs_ij)
}

wij_adjustment<-function(wij,w_t,c,freq_i,freq_j)
{
  w_ij_est<-w_t*(1 - (1 - freq_i)**c)*(1 - (1 - freq_j)**c)
  w_ij_est<-wij-w_ij_est
  w_ij_est[w_ij_est<0]<-0
  w_ij_est
}

derivative_prob_function<-function(f,alpha,W,w,c)
{
  if (f==0){-(alpha*(1-f))+sum(f*c*(W-w)-w)}
  else if (f==1){sum(f*c*(W-w))}
  else {-(alpha*(1-f))+f*sum(c*((W-w) + w*((1-f)**c)/(((1-f)**max(c,0.001)) - 1)))}
}

#get freq est from the well data. 
find_freq<-function(w,w_tot=96,cpw=1000,a=1)
  uniroot(function(x)derivative_prob_function(x,a,w_tot,w,cpw),lower=0,upper=1,tol = 1e-9)$root

get_input<-function(vec,vec2) #vector form: gets two logical vectors, output wi, wj, wij (well data)
{
  wij=sum((vec+vec2)==2)
  wj=sum((vec-vec2)== -1)
  wi=sum((vec2-vec)== -1)
  c(wi,wj,wij)
}

get_input_mat<-function(mat,mat2) #matrix form: gets two matrices outputs wi, wj, wij matrix (well data)
{
  wij=rowSums(mat&mat2)
  wj=rowSums(mat2&!(mat))
  wi=rowSums(mat&!(mat2))
  cbind(wi,wj,wij)
}

