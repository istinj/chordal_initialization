
function v_idx=chordalMatrixIndex(pose_index, num_poses)
  global pose_dim;
  global landmark_dim;
  global flat_rotation_dimension;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  end
  v_idx=1+(pose_index-1)*flat_rotation_dimension;
end
