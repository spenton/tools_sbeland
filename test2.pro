pro test2_1

    ; change date
    vms0 = sd2vms(4110.0)
    vms1 = sd2vms(4115.0)
    sql ="SELECT * FROM dbo.EventsSummary where startTime >='"+vms0+"' and stopTime <='"+vms1+"' and eventTypeName='SolarPeriod' and missionName='SORCE'"
    query_database, sql, solarPeriod, nrows
    help,solarPeriod

end
;
pro test2

    vms0 = sd2vms(4100.0)
    vms1 = sd2vms(4102.0)
    sql ="SELECT * FROM dbo.EventsSummary where startTime >='"+vms0+"' and stopTime <='"+vms1+"' and eventTypeName='SolarPeriod' and missionName='SORCE'"
    query_database, sql, solarPeriod, nrows, server='PS-DB', database='PS', user='anon', pass='anon$password', /dbconnect
    help,solarPeriod

    test2_1

    query_database,/dbclose
end
