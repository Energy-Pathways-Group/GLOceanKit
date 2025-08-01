function mustBeDoulbyPeriodicFPlane(a)
assert( isa(a,'WVGeometryDoublyPeriodic') && isa(a,'WVRotatingFPlane'),'mustBeDoulbyPeriodicFPlane::invalidClass','This geostrophic component is only valid for double periodic geometry on a rotating f-plane.');
end