open Ocamlbuild_plugin;;
open Command;;

dispatch begin function
  | Before_rules -> 
    begin
    end
  | After_rules ->
    begin
      flag ["ocaml";"compile";"native";"gprof"] (S [ A "-p"]);
    end
  | _ -> ()
end
