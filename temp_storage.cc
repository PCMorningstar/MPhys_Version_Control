// ============================================================
// Flavour truth b and bbar indices
// ============================================================
RVec<int> flavour_truth_b_and_bbar_indices(
  const RVec<int>& flavour_indices
){
  return {flavour_indices[4], flavour_indices[6]};
}
// ============================================================
// Flavour truth c index
// ============================================================
int flavour_truth_c_index(
  const RVec<int>& flavour_indices
){
  return flavour_indices[3];
}
// ============================================================
// Flavour light jet indices
// ============================================================
RVec<int> flavour_light_jet_indices(
  const RVec<int>& flavour_indices
){
  return {flavour_indices[0], flavour_indices[1], flavour_indices[2]};
}