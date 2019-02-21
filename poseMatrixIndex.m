
% retrieves the index in the perturbation vector, that corresponds to
% a certain pose
% input:
%   pose_index:     the index of the pose for which we want to compute the
%                   index
%   num_poses:      number of pose variables in the state
%   num_landmarks:  number of pose variables in the state
% output:
%   v_idx: the index of the sub-vector corrsponding to
%          pose_index, in the array of perturbations  (-1 if error)


function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  end
  v_idx=1+(pose_index-1)*pose_dim;
end
