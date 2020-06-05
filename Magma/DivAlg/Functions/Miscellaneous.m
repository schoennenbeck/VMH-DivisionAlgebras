intrinsic RealToString(r::FldReElt) -> MonStgElt
 {This converts a floating point real number into a string.}
 if Sign(r) eq -1 then
  str := "-";
 else
  str := "";
 end if;
 r:=Abs(r);
 p := Integers()! Floor(r) ;
 str := str cat IntegerToString(p) cat ".";
 for i := 1 to 15 do
  r:=10*(r-p);
  p := Integers()! Floor(r) ;
  str := str cat IntegerToString(p);
 end for;
 return str;
end intrinsic;
