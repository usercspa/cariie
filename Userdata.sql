create database Userdata;
use Userdata;


CREATE TABLE IF NOT EXISTS Usersigup(
    
    phone_number VARCHAR(15) NOT NULL,
    country VARCHAR(50) NOT NULL,
    fullname VARCHAR(50) NOT NULL,
    age INT NOT NULL,
    in_period BIT NOT NULL,
    startdate VARCHAR(50) NOT NULL,
    enddate VARCHAR(50) NOT NULL,
    avgperiod_duration INT NOT NULL,
    pain_level VARCHAR(10) NOT NULL,
    bmi FLOAT DEFAULT 0,
    bbt FLOAT DEFAULT 36,
    PRIMARY KEY (phone_number)
);

INSERT INTO Usersigup(phone_number,country,fullname, age, in_period, startdate, enddate, avgperiod_duration,pain_level, bmi, bbt) VALUES
("+12345678902", "Malawi", "James Cameron", 15, 0, '2015-12-10', '2015-12-17', 7, 6, 24, 36.7),
("+12345638902", "Malawi",  "Kim Cameron", 15, 0, '2015-12-10', '2015-12-17', 7, 6, 24, 36.7),
;


CREATE TABLE IF NOT EXISTS Usercollect(
    
    phone_number VARCHAR(15) NOT NULL,
    in_period BIT NOT NULL,
    startdate VARCHAR(50) NOT NULL,
    enddate VARCHAR(50) NOT NULL,
    period_duration INT NOT NULL,
    pain_level VARCHAR(10) NOT NULL,
    bmi FLOAT DEFAULT 0,
    bbt FLOAT DEFAULT 36,
    PRIMARY KEY (phone_number)
);

INSERT INTO Usercollect(phone_number, in_period, startdate, enddate, period_duration, pain_level, bmi, bbt) VALUES
("+12345678902",1, '2016-1-10', '2016-1-17', 7, 6, 24, 36.7)
("+12345638902",1, '2016-1-10', '2016-1-17', 7, 6, 24, 36.7)
;



Create Table Userdata as
select *
from Usersigup, Usercollect
where Usersigup.phone_number = Usercollect.phone_number;


CREATE TABLE Userprofile as
select *
from Userdata
where Userdata.phone_number = "+12345678902";



CREATE TABLE Platform as
select startdate, enddate, period_duration,pain_level, bmi, bbt
from Userdata;