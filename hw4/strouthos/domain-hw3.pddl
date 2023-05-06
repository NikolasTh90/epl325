(define (domain hw3)
    (:requirements :strips :typing :negative-preconditions :disjunctive-preconditions :continuous-effects :equality)
    (:types

        train-station - location
        harbor - location
        warehouse - location
        truck - vehicle
        ship - vehicle
        wagon - vehicle
        train_machine package city country

    )
    (:predicates
        (at ?o - object ?l - location)
        (in_city ?l - location ?c - city)
        (in_country ?c - city ?cn - country)
        (in_vehicle ?p - package ?v - vehicle)
        (AvailableTruckWagon ?v - vehicle)
        (checkShip ?v vehicle)
        (connected_to_train ?w - wagon ?tm - train_machine)
    )

    (:action move-truck
        :parameters (?t - truck ?from - location ?to - location)
        :precondition (and
            (at ?t ?from) ; check if truck at departure
            (forall (?c1 - city)
                (forall (?c2 - city)
                    (imply (and
                            (in_city ?l1 ?c1)
                            (in_city ?l2 ?c2)
                            )
                        (forall (?cn1 - country)
                            (forall (?cn2 - country)
                                (imply (and
                                    (in_country ?c1 ?cn1)
                                    (in_country ?c2 ?cn2)
                                )
                                ((= ?cn1 ?cn2))
                                )
                            )
                        )
                    )
                )
            )            
        )
        :effect (and
                    (not (at ?t ?from)) ; truck is not to location departure
                    (at ?t ?to) ; truck is to destination   
                    (forall(?p - package) ; update the location of the packages carried by the ship
                        (when (in_vehicle ?p ?s)
                            (and
                                (at ?p ?l2)
                                (not(at ?p ?l1))
                            )
                        )
                    )
                )
    ) 


    (:action move_ship
        :parameters (?s - ship ?l1 - location ?l2 - location)
        :precondition (and
            (at ?s ?l1) ; ship is at departure location/harbor
            (forall (?c1 - city)
                (forall (?c2 - city)
                    (imply (and
                            (in_city ?l1 ?c1)
                            (in_city ?l2 ?c2)
                            )
                        (forall (?cn1 - country)
                            (forall (?cn2 - country)
                                (imply (and
                                    (in_country ?c1 ?cn1)
                                    (in_country ?c2 ?cn2)
                                )
                                (not(= ?cn1 ?cn2))
                                )
                            )
                        )
                    )
                )
            )
        )
        :effect (and
            (not (at ?s ?l1)) ; update the location of the ship
            (at ?s ?l2)
            (forall(?p - package) ; update the location of the packages carried by the ship
                (when (in_vehicle ?p ?s)
                    (and
                        (at ?p ?l2)
                        (not(at ?p ?l1))
                    )
                )
            )
            )
        )    
    
    
    (:action move-train
        :parameters (?t - train_machine ?from - location ?to - location)
        :precondition (and
            (at ?t ?from) ; check if train at departure            
            ; (connected-by-train ?from ?to) ; check destination and departure location are connected via train -> done in init
        )
        :effect (and
                    (not (at ?t ?from)) ; train is not to location departure
                    (at ?t ?to) ; train is to destination
                    (forall (?w - wagon) ; for each wagon connected to the train 
                        (when(connected_to_train ?w ?t)
                            (and
                                (forall(?p - package) ; for each package on the wagon
                                    (when(in_vehicle ?p ?w)
                                        (and
                                            (not (at ?p ?from)) ; packages are not to departure
                                            (at ?p ?to) ; packages are now to destination
                                        )
                                    )
                                )
                                (not (at ?w ?to)) ; wagons are not to departure
                                (at ?w ?to) ; wagons now are to destination
                            )
                        )
                    )
                )
    ) 
        
    (:
        (:action loadTruckWagon
            :parameters (?p - package, ?v - vehicle, ?l - location)
            :precondition (and
                (at ?p ?l) ; Check package is in location same with vehicle
                (at ?v ?l)
                (not (checkShip ?v)) ; check not ship the vehicle, the way wagon,truck load is the same
                (not (in_vehicle ?p ?v)) ; check if package not inside truck/wagon
                (AvailableTruckWagon ?v); check  vehicle truck/wagon is available -> true
            )
            :effect (and
                (in_vehicle ?p ?v)
                (not (AvailableTruckWagon ?v)) ; make  vehicle truck/wagon not available
            )
        )
    ) 

    (:
        (:action unloadTruckWagon
            :parameters (?p - package, ?v - vehicle) ; maybe add location
            :precondition (and
                (in_vehicle ?p ?v)
                (not (AvailableTruckWagon ?v)) ; check  vehicle truck/wagon is not available
            )
            :effect (and
                (not (in_vehicle ?p ?v)) ; package not inside truck/wagon
                (AvailableTruckWagon ?v) ; check  vehicle truck/wagon is empty
            )
        )
    ) 

    (:
        (:action loadShip
            :parameters (?p - package, ?v - vehicle, ?l - location)
            :precondition (and
                (at ?p ?l) ; Check package is in location same with vehicle
                (at ?v ?l)
                ((checkShip ?v)) ; check if ship
                (not (in_vehicle ?p ?v)) ; check if package not inside ship
            )
            :effect (and
                (in_vehicle ?p ?v) ; package is in the ship
            )
        )
    )
     
    (:
        (:action unloadShip
            :parameters (?p - package, ?v - vehicle) 
            :precondition (and
                (in_vehicle ?p ?v) ; package inside ship
            )
            :effect (and
                (not (in_vehicle ?p ?v)) ; package not inside ship 
            )
        )
    ) 

    (:action connect-wagon
        :parameters (?w - wagon ?tm - train_machine ?l - location)
        :precondition (and 
            (at ?tm ?l) ; check that wagon and tm is to the same location 
            (at ?w ?l)
            (forall (?tm - train_machine)  ; check that wagon is not connected to any other train_machine
                (and
                    ((not)connected_to_train ?w ?tm)
                )
            ) 
        ) 
        :effect (and             
            (connected_to_train ?w ?tm)

        )
    )
        
    (:action unconnect-wagon
        :parameters (?w - wagon ?tm - train_machine ?l - location)
        :precondition (and 
            (at ?tm ?l) ; check that wagon and tm is to the same location 
            (at ?w ?l)
            (connected_to_train ?w ?tm)
        ) 
        :effect (and             
            (not (connected_to_train ?w ?tm))

        )
    )

)