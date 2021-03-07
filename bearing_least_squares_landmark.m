1;

function [e,Jr, Jl]=landmarkErrorAndJacobian(Xr,Xl,z)
  R = Xr(1:2,1:2);
  t = Xr(1:2,3);
  p_hat = R'*(Xl-t);
  z_hat = atan2(p_hat(2),p_hat(1));
  e = z_hat-z;
  e = atan2(sin(e), cos(e)); 
  J_r = zeros(2,3);
  J_l = zeros(2,2);
  dR0_T = [0 1;-1 0];
  J_r(1:2,1:2) = -R';               
  J_r(1:2,3) = R'*dR0_T*Xl;     
  J_l(:,:) = R';                    
  J_pre = (1./(p_hat(1)^2+p_hat(2)^2))*[-p_hat(2) p_hat(1)];           
  Jr = zeros(1,3);
  Jl = zeros(1,2);
  Jr = J_pre*J_r;                  
  Jl = J_pre*J_l;                    
endfunction;

function [H, b, chi_tot, num_inliers]=buildLinearSystemLandmarks(XR,XL,ZL,landmark_associations,kernel_threshold)
  global pose_dim;
  global landmark_dim;
  num_poses = size(XR,3);
  num_landmarks = size(XL,2);

  system_size = pose_dim*num_poses+landmark_dim*num_landmarks; 
  H=zeros(system_size,system_size);
  b=zeros(system_size,1);
  chi_tot=0;
  chi=0;
  num_inliers=0;
  for i=1:size(ZL,2)
    pose_index = landmark_associations(1,i);
    landmark_index = landmark_associations(2,i);
    Xr = XR(:,:,pose_index);
    Xl = XL(:,landmark_index);
    z = ZL(1,i);
    [e,Jr, Jl]=landmarkErrorAndJacobian(Xr,Xl,z);
    chi=e'*e;
    if (chi>kernel_threshold)
      e*=sqrt(kernel_threshold/chi);
      chi=kernel_threshold;
    else
      num_inliers++;
    endif;
    chi_tot+=chi;
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

    H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jr'*Jr;

    H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jr'*Jl;

    H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jl'*Jl;

    H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jl'*Jr;

    b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jr'*e;
      b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jl'*e;
    endfor;
endfunction;