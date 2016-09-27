 function [ c r ] = randomSphere(dims, radius)
 % creating one sphere at random inside [0..dims(1)]x[0..dims(2)]x...
 % radius and center coordinates are sampled from a uniform distribution 
 % over the relevant domain.

 r = random('normal',radius(1),radius(2));
 c = bsxfun(@times,(dims - 2*r) , rand(1,3)) + r; 
 end