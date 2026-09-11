# ============================================================
# 06_prepare_weather_calendar.R
# Melbourne Pedestrian Bayesian Project
#
# Prepare 2025 calendar variables and Victorian public holidays
# ============================================================

library(here)
library(dplyr)
library(lubridate)
library(tibble)

calendar_2025 <- tibble(
  date = seq.Date(
    from = as.Date("2025-01-01"),
    to   = as.Date("2025-12-31"),
    by   = "day"
  )
) %>%
  mutate(
    day_of_week = wday(
      date,
      label = TRUE,
      abbr = FALSE,
      week_start = 1
    ),
    
    weekend = wday(
      date,
      week_start = 1
    ) >= 6,
    
    month = month(
      date,
      label = TRUE,
      abbr = FALSE
    ),
    
    season = case_when(
      month(date) %in% c(12, 1, 2) ~ "Summer",
      month(date) %in% c(3, 4, 5)  ~ "Autumn",
      month(date) %in% c(6, 7, 8)  ~ "Winter",
      month(date) %in% c(9, 10, 11) ~ "Spring"
    )
  )

vic_public_holidays_2025 <- tibble(
  date = as.Date(c(
    "2025-01-01",
    "2025-01-27",
    "2025-03-10",
    "2025-04-18",
    "2025-04-19",
    "2025-04-20",
    "2025-04-21",
    "2025-04-25",
    "2025-06-09",
    "2025-09-26",
    "2025-11-04",
    "2025-12-25",
    "2025-12-26"
  )),
  
  holiday_name = c(
    "New Year's Day",
    "Australia Day",
    "Labour Day",
    "Good Friday",
    "Saturday before Easter Sunday",
    "Easter Sunday",
    "Easter Monday",
    "ANZAC Day",
    "King's Birthday",
    "Friday before AFL Grand Final",
    "Melbourne Cup Day",
    "Christmas Day",
    "Boxing Day"
  )
)

calendar_2025 <- calendar_2025 %>%
  left_join(
    vic_public_holidays_2025,
    by = "date"
  ) %>%
  mutate(
    public_holiday = !is.na(holiday_name)
  )



