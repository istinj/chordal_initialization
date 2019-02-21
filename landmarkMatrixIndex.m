
% retrieves the index in the perturbation vector, that corresponds to
% a certain landmark
% input:
%   landmark_index:     the index of the landmark for which we want to compute the
%                   index
%   num_poses:      number of pose variables in the state
%   num_landmarks:  number of pose variables in the state
% output:
%   v_idx: the index of the perturnation corrsponding to the
%           landmark_index, in the array of perturbations

function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  end
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
end