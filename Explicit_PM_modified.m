function diff_im = Explicit_PM_modified(im, num_iter, delta_t, kappa, option, sigma)
%       ARGUMENT DESCRIPTION:
%               im       - gray scale image (MxN).
%               num_iter - number of iterations. 
%               delta_t  - integration constant (0 <= delta_t <= 0.25).

%               kappa    - gradient modulus threshold that controls the conduction.
%               option  - conduction coefficient functions proposed by Perona & Malik:
%                          1 - c(x,y,t) = exp(-(nablaI/kappa).^2),
%                              privileges high-contrast edges over low-contrast ones. 
%                          2 - c(x,y,t) = 1./(1 + (nablaI/kappa).^2),
%                              privileges wide regions over smaller ones. 
%               sigma   - variance of the gaussian function convolved with
%                         the image to find the values of c(x,y,t)
%       OUTPUT DESCRIPTION:
%                diff_im - (diffused) image with the largest scale-space parameter.

% Convert input image to double.
im = double(im);

% PDE (partial differential equation) initial condition.
diff_im = im;

% Center pixel distances.
dx = 1;
dy = 1;

% 2D convolution masks - finite differences.
 hN = [0 1 0; 0 -1 0; 0 0 0];
 hS = [0 0 0; 0 -1 0; 0 1 0];
 hE = [0 0 0; 0 -1 1; 0 0 0];
 hW = [0 0 0; 1 -1 0; 0 0 0];
%///////////////////////////////////////////////
% hNE = [0 0 1; 0 -1 0; 0 0 0];
% hSE = [0 0 0; 0 -1 0; 0 0 1];
% hSW = [0 0 0; 0 -1 0; 1 0 0];
% hNW = [1 0 0; 0 -1 0; 0 0 0];

% Anisotropic diffusion.
for t = 1:num_iter
   
    % Finite differencies in 4 different directions
    nablaN = imfilter(diff_im,hN,'conv');
    nablaS = imfilter(diff_im,hS,'conv');   
    nablaW = imfilter(diff_im,hE,'conv');
    nablaE = imfilter(diff_im,hW,'conv');
%     nablaNE = imfilter(diff_im,hNE,'conv');
%     nablaSE = imfilter(diff_im,hSE,'conv');   
%     nablaSW = imfilter(diff_im,hSW,'conv');
%     nablaNW = imfilter(diff_im,hNW,'conv');
    % Neumann boundary condition: isolated system
    nablaN(1,:) = 0;
    nablaS(end,:) = 0;
    nablaW(:,1) = 0;
    nablaE(:,end) = 0;
    
%     nablaNE(1,:) = 0;
%     nablaSE(end,:) = 0;
%     nablaSW(:,1) = 0;
%     nablaNW(:,end) = 0;
    if sigma > 0  
        
        diff_blur = imgaussfilt(diff_im,sigma);
        nablaN_blur = imfilter(diff_blur,hN,'conv');
        nablaS_blur = imfilter(diff_blur,hS,'conv');   
        nablaW_blur = imfilter(diff_blur,hE,'conv');
        nablaE_blur = imfilter(diff_blur,hW,'conv');
        
%         nablaNE_blur = imfilter(diff_blur,hNE,'conv');
%         nablaSE_blur = imfilter(diff_blur,hSE,'conv');   
%         nablaSW_blur = imfilter(diff_blur,hSW,'conv');
%         nablaNW_blur = imfilter(diff_blur,hNW,'conv');
        % Neumann boundary condition: isolated system
        nablaN_blur(1,:) = 0;
        nablaS_blur(end,:) = 0;
        nablaW_blur(:,1) = 0;
        nablaE_blur(:,end) = 0; 
%         nablaNE_blur(1,:) = 0;
%         nablaSE_blur(end,:) = 0;
%         nablaSW_blur(:,1) = 0;
%         nablaNW_blur(:,end) = 0;
        
        if option == 1
            cN = exp(-(nablaN_blur/kappa).^2);
            cS = exp(-(nablaS_blur/kappa).^2);
            cW = exp(-(nablaW_blur/kappa).^2);
            cE = exp(-(nablaE_blur/kappa).^2);
%             cNE = exp(-(nablaNE_blur/kappa).^2);
%             cSE = exp(-(nablaSE_blur/kappa).^2);
%             cSW = exp(-(nablaSW_blur/kappa).^2);
%             cNW = exp(-(nablaNW_blur/kappa).^2);
        
        elseif option == 2
            cN = 1./(1 + (nablaN_blur/kappa).^2);
            cS = 1./(1 + (nablaS_blur/kappa).^2);
            cW = 1./(1 + (nablaW_blur/kappa).^2);
            cE = 1./(1 + (nablaE_blur/kappa).^2);
%             
%             cNE = 1./(1 + (nablaNE_blur/kappa).^2);
%             cSE = 1./(1 + (nablaSE_blur/kappa).^2);
%             cSW = 1./(1 + (nablaSW_blur/kappa).^2);
%             cNW = 1./(1 + (nablaNW_blur/kappa).^2);
%             
        
        elseif option==3
            cN=1-exp(-3.315./((nablaN_blur).^4/kappa.^8));
            cS=1-exp(-3.315./((nablaS_blur).^4/kappa.^8));
            cW=1-exp(-3.315./((nablaW_blur).^4/kappa.^8));
            cE=1-exp(-3.315./((nablaE_blur).^4/kappa.^8));
%             cNE=1-exp(-3.315./((nablaNE_blur).^4/kappa.^8));
%             cSE=1-exp(-3.315./((nablaSE_blur).^4/kappa.^8));
%             cSW=1-exp(-3.315./((nablaSW_blur).^4/kappa.^8));
%             cNW=1-exp(-3.315./((nablaNW_blur).^4/kappa.^8));
        end   
    else
        if option == 1
            cN = exp(-(nablaN/kappa).^2);
            cS = exp(-(nablaS/kappa).^2);
            cW = exp(-(nablaW/kappa).^2);
            cE = exp(-(nablaE/kappa).^2);
%             cNE = exp(-(nablaNE/kappa).^2);
%             cSE = exp(-(nablaSE/kappa).^2);
%             cSW = exp(-(nablaSW/kappa).^2);
%             cNW = exp(-(nablaNW/kappa).^2);
            
        elseif option == 2
            cN = 1./(1 + (nablaN/kappa).^2);
            cS = 1./(1 + (nablaS/kappa).^2);
            cW = 1./(1 + (nablaW/kappa).^2);
            cE = 1./(1 + (nablaE/kappa).^2);
%             cNE = 1./(1 + (nablaNE/kappa).^2);
%             cSE = 1./(1 + (nablaSE/kappa).^2);
%             cSW = 1./(1 + (nablaSW/kappa).^2);
%             cNW = 1./(1 + (nablaNW/kappa).^2);  
         
        elseif option==3
           cN=1-exp(-3.315./((nablaN).^4/kappa.^8)); 
           cS=1-exp(-3.315./((nablaS).^4/kappa.^8));
           cW=1-exp(-3.315./((nablaW).^4/kappa.^8));
           cE=1-exp(-3.315./((nablaE).^4/kappa.^8));
%            cNE=1-exp(-3.315./((nablaNE).^4/kappa.^8)); 
%            cSE=1-exp(-3.315./((nablaSE).^4/kappa.^8));
%            cSW=1-exp(-3.315./((nablaSW).^4/kappa.^8));
%            cNW=1-exp(-3.315./((nablaNW).^4/kappa.^8));
        end
    
    end     
    % Discrete PDE solution.
    diff_im = diff_im + ...
        delta_t*(...
        (1/(dy^2))*cN.*nablaN + (1/(dy^2))*cS.*nablaS + ...
        (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE); 
%         (1/(dy^2))*cNE.*nablaNE + (1/(dy^2))*cSW.*nablaSW+...
%         (1/(dx^2))*cNW.*nablaNW + (1/(dx^2))*cSE.*nablaSE
    
    % Iteration warning.
    fprintf('\rIteration %d\n',t);
end

%figure('Name','Explicit Perona Malik');
%imshow(uint8(diff_im));

