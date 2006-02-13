// Copyright Vladimir Prus 2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_PARSERS_HPP_VP_2004_05_06
#define BOOST_PARSERS_HPP_VP_2004_05_06

#include <boost/program_options/detail/convert.hpp>

namespace boost { namespace program_options {

    namespace detail {
        template<class charT, class Iterator>
        std::vector<std::basic_string<charT> > 
        make_vector(Iterator i, Iterator e)
        {
            std::vector<std::basic_string<charT> > result;
            // Some compilers don't have templated constructor for 
            // vector, so we can't create vector from (argv+1, argv+argc) range
            for(; i != e; ++i)
                result.push_back(*i);
            return result;            
        }
    }

    template<class charT>
    basic_command_line_parser<charT>::
    basic_command_line_parser(const std::vector<
                              std::basic_string<charT> >& args)
    : common_command_line_parser(to_internal(args))
    {}


    template<class charT>
    basic_command_line_parser<charT>::
    basic_command_line_parser(int argc, charT* argv[])
    : common_command_line_parser(
        // Explicit template arguments are required by gcc 3.3.1 
        // (at least mingw version), and do no harm on other compilers.
        to_internal(detail::make_vector<charT, charT**>(argv+1, argv+argc)))
    {}

    
    template<class charT>
    basic_command_line_parser<charT>& 
    basic_command_line_parser<charT>::options(const options_description& desc)
    {
        m_desc = &desc;
        return *this;
    }

    template<class charT>
    basic_command_line_parser<charT>& 
    basic_command_line_parser<charT>::positional(
        const positional_options_description& desc)
    {
        m_positional = &desc;
        return *this;
    }

    template<class charT>
    basic_command_line_parser<charT>& 
    basic_command_line_parser<charT>::style(int style)
    {
        m_style = style;
        return *this;
    }

    template<class charT>
    basic_command_line_parser<charT>& 
    basic_command_line_parser<charT>::extra_parser(ext_parser ext)
    {
        m_ext = ext;
        return *this;
    }

    template<class charT>    
    basic_parsed_options<charT>
    basic_command_line_parser<charT>::run() const
    {
        // Presense of parsed_options -> wparsed_options conversion
        // does the trick.
        return basic_parsed_options<charT>(
            common_command_line_parser::run());
    }


    template<class charT>
    basic_parsed_options<charT>
    parse_command_line(int argc, charT* argv[],
                       const options_description& desc,
                       int style,
                       function1<std::pair<std::string, std::string>, 
                                 const std::string&> ext)
    {
        return basic_command_line_parser<charT>(argc, argv).options(desc).
            style(style).extra_parser(ext).run();
    }

}}

#endif
