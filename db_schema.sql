-- phpMyAdmin SQL Dump
-- version 3.4.10.1deb1
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Sep 22, 2014 at 02:44 PM
-- Server version: 5.5.37
-- PHP Version: 5.3.10-1ubuntu3.11

SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;

--
--

-- --------------------------------------------------------

--
-- Table structure for table `blast6out`
--

CREATE TABLE IF NOT EXISTS `blast6out` (
  `id` int(11) NOT NULL,
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL,
  PRIMARY KEY (`id`),
  KEY `evalue` (`evalue`),
  KEY `percent_identity` (`percent_identity`),
  KEY `alignment_length` (`alignment_length`),
  KEY `bitscore` (`bitscore`),
  FULLTEXT KEY `query_label` (`query_label`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lessstrict_0`
--

CREATE TABLE IF NOT EXISTS `lessstrict_0` (
  `id` int(11) NOT NULL DEFAULT '0',
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lessstrict_1`
--

CREATE TABLE IF NOT EXISTS `lessstrict_1` (
  `id` int(11) NOT NULL DEFAULT '0',
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lessstrict_2`
--

CREATE TABLE IF NOT EXISTS `lessstrict_2` (
  `id` int(11) NOT NULL DEFAULT '0',
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lessstrict_3`
--

CREATE TABLE IF NOT EXISTS `lessstrict_3` (
  `id` int(11) NOT NULL DEFAULT '0',
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lessstrict_4`
--

CREATE TABLE IF NOT EXISTS `lessstrict_4` (
  `id` int(11) NOT NULL DEFAULT '0',
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lessstrict_5`
--

CREATE TABLE IF NOT EXISTS `lessstrict_5` (
  `id` int(11) NOT NULL DEFAULT '0',
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Table structure for table `lessstrict_6`
--

CREATE TABLE IF NOT EXISTS `lessstrict_6` (
  `id` int(11) NOT NULL DEFAULT '0',
  `query_label` varchar(255) NOT NULL,
  `location` int(2) DEFAULT NULL,
  `target` varchar(255) NOT NULL,
  `percent_identity` double NOT NULL,
  `alignment_length` int(11) NOT NULL,
  `mismatch` int(11) NOT NULL,
  `gap` int(11) NOT NULL,
  `query_start` int(11) NOT NULL,
  `query_end` int(11) NOT NULL,
  `target_start` int(11) NOT NULL,
  `target_end` int(11) NOT NULL,
  `evalue` double NOT NULL,
  `bitscore` double NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;

