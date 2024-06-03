/**
 * @file demand.h, part of the project TransOMS under Apache License 2.0
 * @author jdlph (jdlph@hotmail.com) and xzhou99 (xzhou74@asu.edu)
 * @brief Definitions of classes related to demand
 *
 * @copyright Copyright (c) 2023 Peiheng Li, Ph.D. and Xuesong (Simon) Zhou, Ph.D.
 */

#ifndef GUARD_DEMAND_H
#define GUARD_DEMAND_H

#include <global.h>

#include <map>
#include <memory>
#include <vector>

namespace transoms
{
class Agent {
public:
    Agent() = delete;

    Agent(size_type no_, unsigned short at_no_, unsigned short dp_no_,
          size_type oz_no_, size_type dz_no_, const Column* c = nullptr)
        : no {no_}, at_no {at_no_}, dp_no {dp_no_}, oz_no {oz_no_}, dz_no {dz_no_},
          col {c}, pce {1}
    {
        initialize_intervals();
    }

    Agent(const Agent&) = delete;
    Agent& operator=(const Agent&) = delete;

    Agent(Agent&&) = default;
    Agent& operator=(Agent&&) = delete;

    ~Agent() = default;

    bool completes_trip() const
    {
        return dep_intvls.front() > 0;
    }

    auto get_agent_type_no() const
    {
        return at_no;
    }

    const Column* get_column() const
    {
        return col;
    }

    auto get_demand_period_no() const
    {
        return dp_no;
    }

    auto get_dest_zone_no() const
    {
        return dz_no;
    }

    auto get_orig_zone_no() const
    {
        return oz_no;
    }

    auto get_od() const
    {
        return std::make_pair(oz_no, dz_no);
    }

    size_type get_no() const
    {
        return no;
    }

    double get_pce() const
    {
        return pce;
    }

    // simulation
    bool reaches_last_link() const
    {
        return curr_link_no == 0;
    }

    auto get_arr_interval() const
    {
        return arr_intvls[curr_link_no];
    }

    auto get_dep_interval() const
    {
        return dep_intvls[curr_link_no];
    }

    size_type get_next_link_no() const
    {
        return get_link_path()[curr_link_no - 1];
    }

    size_type get_dest_arr_interval() const
    {
        auto i = get_dep_interval();
        if (i < std::numeric_limits<size_type>::max())
            return i;

        return dep_intvls[curr_link_no + 1];
    }

    double get_orig_dep_time() const
    {
        return dep_time;
    }

    size_type get_orig_dep_interval() const
    {
        return dep_intvls.back();
    }

    size_type get_travel_interval() const
    {
        return get_dest_arr_interval() - arr_intvls.back();
    }

    void move_to_next_link()
    {
        if (curr_link_no > 0)
            --curr_link_no;
    }

    void set_arr_interval(unsigned i, size_type increment = 0)
    {
        arr_intvls[curr_link_no - increment] = i;
    }

    void set_dep_interval(unsigned i)
    {
        dep_intvls[curr_link_no] = i;
    }

    void set_dep_time(double t)
    {
        dep_time = t;
    }

    // can be combined with set_arr_interval()?
    void set_orig_arr_interval(unsigned i)
    {
        arr_intvls.back() = i;
    }

    void increment_dep_interval(size_type i)
    {
        dep_intvls[curr_link_no] = arr_intvls[curr_link_no] + i;
    }

    const std::vector<size_type>& get_link_path() const;

    std::vector<size_type>::size_type get_link_num() const;
    size_type get_first_link_no() const;

    std::vector<size_type> get_time_sequence() const;

private:
    void initialize_intervals();

private:
    size_type no;

    unsigned short at_no;
    unsigned short dp_no;

    // use unsigned short instead?
    size_type oz_no;
    size_type dz_no;

    const Column* col;

    double pce;

    // simulation
    size_type curr_link_no;
    double dep_time;

    std::vector<size_type> arr_intvls;
    std::vector<size_type> dep_intvls;
};

class AgentType {
public:
    AgentType() : no {0}, name {"auto"}, flow_type {0}, pce {1},
                  vot {10}, ffs {60}, is_link_ffs {true}
    {
    }

    AgentType(unsigned short no_, std::string& name_, unsigned short flow_type_,
              double pce_, double vot_,  double ffs_, bool use_link_ffs_)
        : no {no_}, name {std::move(name_)}, flow_type {flow_type_}, pce {pce_},
          vot {vot_}, ffs {ffs_}, is_link_ffs {use_link_ffs_}
    {
    }

    AgentType(const AgentType&) = delete;
    AgentType& operator=(const AgentType&) = delete;

    AgentType(AgentType&&) = delete;
    AgentType& operator=(AgentType&&) = delete;

    ~AgentType() = default;

    auto get_no() const
    {
        return no;
    }

    auto get_flow_type() const
    {
        return flow_type;
    }

    auto get_ffs() const
    {
        return ffs;
    }

    const std::string& get_name() const
    {
        return name;
    }

    auto get_pce() const
    {
        return pce;
    }

    auto get_vot() const
    {
        return vot;
    }

    bool use_link_ffs() const
    {
        return is_link_ffs;
    }

public:
    static const std::string& get_default_name()
    {
        return AT_DEFAULT_NAME;
    }

    static const std::string& get_legacy_name()
    {
        return AT_LEGACY_NAME;
    }

private:
    unsigned short no;
    std::string name;

    unsigned short flow_type;
    double pce;
    double vot;

    double ffs;
    bool is_link_ffs;
};

class Demand {
public:
    Demand() = delete;

    explicit Demand(const AgentType* at_) : at {at_}
    {
    }

    Demand(unsigned short no_, std::string& filename_, const AgentType* at_)
        : no {no_}, filename {std::move(filename_)}, at {at_}
    {
    }

    Demand(const Demand&) = default;
    Demand& operator=(const Demand&) = delete;

    Demand(Demand&&) = default;
    Demand& operator=(Demand&&) = delete;

    ~Demand() = default;

    auto get_no() const
    {
        return no;
    }

    const std::string& get_agent_type_name() const
    {
        return at->get_name();
    }

    auto get_agent_type_no() const
    {
        return at->get_no();
    }

    const std::string& get_file_name() const
    {
        return filename;
    }

private:
    unsigned short no = 0;
    std::string filename = "demand.csv";

    const AgentType* at;
};

class DemandPeriod {
public:
    DemandPeriod() : no {0}, period {"AM"}, time_period {"0700-0800"}, se {nullptr}
    {
    }

    explicit DemandPeriod(Demand&& dem) : DemandPeriod()
    {
        ds.push_back(dem);
    }

    DemandPeriod(unsigned short no_,
                 std::string& period_, std::string& time_period_,
                 Demand&& dem, std::unique_ptr<SpecialEvent>& se_)
        : no {no_}, period {std::move(period_)}, time_period {std::move(time_period_)}, se {std::move(se_)}
    {
        ds.push_back(dem);
        setup_time();
    }

    DemandPeriod(const DemandPeriod&) = delete;
    DemandPeriod& operator=(const DemandPeriod&) = delete;

    DemandPeriod(DemandPeriod&&) = delete;
    DemandPeriod& operator=(DemandPeriod&&) = delete;

    ~DemandPeriod() = default;

    auto get_no() const
    {
        return no;
    }

    const std::string& get_period() const
    {
        return period;
    }

    const std::string& get_time_period() const
    {
        return time_period;
    }

    const auto& get_demands() const
    {
        return ds;
    }

    const auto& get_special_event() const
    {
        return se;
    }

    // minute as time of day
    unsigned short get_start_time() const
    {
        return start_time;
    }

    unsigned short get_end_time() const
    {
        return start_time + dur;
    }

    unsigned short get_duration() const
    {
        return dur;
    }

    bool contain_iter_no(unsigned short iter_no) const;
    double get_cap_ratio(const std::string& link_id, unsigned short iter_no) const;

private:
    void setup_time();
    unsigned short to_minutes(const std::string& t);

private:
    unsigned short no;

    std::string period;
    std::string time_period;

    unsigned short start_time = 420;
    unsigned short dur = 60;

    std::vector<Demand> ds;
    const std::unique_ptr<SpecialEvent> se;
};

} // namespace transoms

#endif