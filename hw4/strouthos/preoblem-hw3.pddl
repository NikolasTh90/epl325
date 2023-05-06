(define (problem problem-hw3-ex1)
    (:domain hw3)
    (:objects
        Spain Italy - country
        Barcelona Madrid Rome - city
        BarcaHB MadHB RomeHB - harbor
        MadridST RomeST - train-station
        BarcaWH BarcaWH2 MadWH RomeWH - warehouse
        ship1 ship2 - ship
        truck1 truck2 - truck
        tm1 tm2 - train_machine
        wagon1 wagon2 wagon3 wagon4 wagon5 wagon6 - wagon
        package1 package2 package3 package4 package5 package6 package7 package8 package9 package10 - package
    )
    (:init
        (at ship1 BarcaHB)
        (at ship2 MadHB)
        (at truck1 BarcaWH)
        (at truck2 MadWH)
        (at tm1 MadridST)
        (at tm2 RomeST)
        (at wagon1 MadridST)
        (at wagon2 MadridST)
        (at wagon3 RomeST)
        (at wagon4 RomeST)
        (at wagon5 BarcaHB)
        (at wagon6 BarcaHB)
        (at package1 BarcaHB)
        (at package2 MadHB)
        (at package3 MadHB)
        (at package4 MadridST)
        (at package5 RomeST)
        (at package6 BarcaWH)
        (at package7 BarcaWH)
        (at package8 BarcaWH2)
        (at package9 MadWH)
        (at package10 RomeWH)
        (connected_to_train wagon1 tm1)
        (connected_to_train wagon2 tm1)
        (connected_to_train wagon3 tm2)
        (connected_to_train wagon4 tm2)

        (in_city BarcaHB Barcelona)
        (in_city BarcaWH Barcelona)
        (in_city BarcaWH2 Barcelona)
        (in_city MadHB Madrid)
        (in_city MadridST Madrid)
        (in_city MadWH Madrid)
        (in_city RomeHB Rome)
        (in_city RomeST Rome)
        (in_city RomeWH Rome)
        (in_country Barcelona Spain)
        (in_country Madrid Spain)
        (in_country Rome Italy)
        (AvailableTruckWagon truck1)
        (AvailableTruckWagon truck1)
        (AvailableTruckWagon wagon1)
        (AvailableTruckWagon wagon2)
        (AvailableTruckWagon wagon3)
        (AvailableTruckWagon wagon4)
        (AvailableTruckWagon wagon5)
        (AvailableTruckWagon wagon6)
        

    )
    (:goal (and
        (at package1 BarcaWH)
        (at package2 RomeWH)
        (at package3 MadWH)
        (at package4 BarcaHB)
        (at package5 MadHB)
        (at package6 RomeST)
        (at package7 RomeST)
        (at package8 RomeST)
        (at package9 RomeST)
        (at package10 RomeST)
    ))
)